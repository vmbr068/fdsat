#ifndef OQPSK_OFFLINE
#include <iio.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace oqpsk {
constexpr double pi = 3.14159265358979323846;
constexpr long long sample_rate = 2'400'000;
constexpr long long symbol_rate = 10'000;
constexpr double rolloff = 0.25;
constexpr int sps = static_cast<int>(sample_rate / symbol_rate);
constexpr int span = 16;
[[maybe_unused]] constexpr long long tx_frequency = 433'000'000;
[[maybe_unused]] constexpr long long rf_bandwidth = 200'000;
[[maybe_unused]] constexpr double tx_gain_db = -30.0;
constexpr size_t buffer_samples = 65'536;
constexpr float peak_limit = 0.8f;
static_assert(sample_rate % symbol_rate == 0 && sps % 2 == 0);
struct IQ {
  float i;
  float q;
};
static_assert(sizeof(IQ) == 2 * sizeof(float) && sizeof(float) == 4);
volatile std::sig_atomic_t stop_requested = 0;
void signal_handler(int) { stop_requested = 1; }

std::vector<float> design_rrc_filter(double beta, int samples_per_symbol,
                                     int span_symbols) {
  if (!(beta > 0.0 && beta <= 1.0) || samples_per_symbol < 2 ||
      samples_per_symbol % 2 || span_symbols < 2 || span_symbols % 2)
    throw std::runtime_error("Invalid RRC parameters");
  const size_t count =
      static_cast<size_t>(span_symbols) * samples_per_symbol + 1;
  std::vector<float> taps(count);
  for (size_t n = 0; n < count; ++n) {
    const double t =
        (static_cast<double>(n) - (count - 1) / 2.0) / samples_per_symbol;
    double h;
    if (std::abs(t) < 1e-12) {
      h = 1.0 + beta * (4.0 / pi - 1.0);
    } else if (std::abs(std::abs(4.0 * beta * t) - 1.0) < 1e-10) {
      h = beta / std::sqrt(2.0) *
          ((1.0 + 2.0 / pi) * std::sin(pi / (4.0 * beta)) +
           (1.0 - 2.0 / pi) * std::cos(pi / (4.0 * beta)));
    } else {
      h = (std::sin(pi * t * (1.0 - beta)) +
           4.0 * beta * t * std::cos(pi * t * (1.0 + beta))) /
          (pi * t * (1.0 - 16.0 * beta * beta * t * t));
    }
    taps[n] = static_cast<float>(h);
  }
  return taps;
}

class RandomPayload {
public:
  explicit RandomPayload(uint32_t seed) : rng_(seed) {}
  IQ next_symbol() {
    if (cursor_ == bits_.size()) {
      for (auto &bit : bits_)
        bit = static_cast<uint8_t>(bit_(rng_));
      cursor_ = 0;
    }
    const float i = bits_[cursor_++] ? 1.0f : -1.0f;
    const float q = bits_[cursor_++] ? 1.0f : -1.0f;
    return {i, q};
  }

private:
  std::mt19937 rng_;
  std::uniform_int_distribution<int> bit_{0, 1};
  std::array<uint8_t, 1024> bits_{};
  size_t cursor_ = bits_.size();
};

class Modulator {
public:
  explicit Modulator(uint32_t seed)
      : payload_(seed), taps_(design_rrc_filter(rolloff, sps, span)) {
    double bound = 0.0;
    for (int phase = 0; phase < sps; ++phase) {
      double sum = 0.0;
      for (size_t k = phase; k < taps_.size(); k += sps)
        sum += std::abs(taps_[k]);
      bound = std::max(bound, sum);
    }
    gain_ = static_cast<float>(peak_limit / (std::sqrt(2.0) * bound));
  }

  IQ next() {
    if (phase_ == 0) {
      for (size_t k = history_.size() - 1; k > 0; --k)
        history_[k] = history_[k - 1];
      history_[0] = payload_.next_symbol();
    }
    IQ shaped{0.0f, 0.0f};
    size_t symbol = 0;
    for (size_t k = phase_; k < taps_.size(); k += sps, ++symbol) {
      shaped.i += history_[symbol].i * taps_[k];
      shaped.q += history_[symbol].q * taps_[k];
    }
    const float delayed_q = q_delay_[q_cursor_];
    q_delay_[q_cursor_] = shaped.q;
    q_cursor_ = (q_cursor_ + 1) % q_delay_.size();
    phase_ = (phase_ + 1) % sps;
    return {gain_ * shaped.i, gain_ * delayed_q};
  }

  void fill(IQ *output, size_t count) {
    for (size_t n = 0; n < count; ++n)
      output[n] = next();
  }
  void prime() {
    for (int n = 0; n < (span + 2) * sps; ++n)
      (void)next();
  }
  float gain() const { return gain_; }

private:
  RandomPayload payload_;
  std::vector<float> taps_;
  std::array<IQ, span + 1> history_{};
  std::array<float, sps / 2> q_delay_{};
  size_t q_cursor_ = 0;
  int phase_ = 0;
  float gain_ = 0.0f;
};

#ifndef OQPSK_OFFLINE
void check_iio_error(long long result, const std::string &operation) {
  if (result < 0) {
    char message[256]{};
    iio_strerror(static_cast<int>(-result), message, sizeof(message));
    throw std::runtime_error(operation + ": " + message);
  }
}

class PlutoTransmitter {
public:
  explicit PlutoTransmitter(const std::string &uri) {
    try {
      ctx_ = iio_create_context_from_uri(uri.c_str());
      if (!ctx_)
        throw std::runtime_error("Failed to create IIO context: " + uri);
      check_iio_error(iio_context_set_timeout(ctx_, 2000), "IIO timeout");
      auto *phy = iio_context_find_device(ctx_, "ad9361-phy");
      auto *tx = iio_context_find_device(ctx_, "cf-ad9361-dds-core-lpc");
      if (!phy || !tx)
        throw std::runtime_error("Pluto devices not found");
      auto *lo = iio_device_find_channel(phy, "altvoltage1", true);
      auto *phy_tx = iio_device_find_channel(phy, "voltage0", true);
      if (!lo || !phy_tx)
        throw std::runtime_error("PHY TX channels not found");
      check_iio_error(iio_channel_attr_write(phy_tx, "rf_port_select", "A"),
                      "TX port");
      check_iio_error(iio_channel_attr_write_longlong(
                          phy_tx, "sampling_frequency", sample_rate),
                      "TX sample rate");
      check_iio_error(
          iio_channel_attr_write_longlong(phy_tx, "rf_bandwidth", rf_bandwidth),
          "RF bandwidth");
      check_iio_error(
          iio_channel_attr_write_double(phy_tx, "hardwaregain", tx_gain_db),
          "TX gain");
      check_iio_error(
          iio_channel_attr_write_longlong(lo, "frequency", tx_frequency),
          "TX frequency");
      long long actual_fs = 0, actual_bw = 0;
      check_iio_error(iio_channel_attr_read_longlong(
                          phy_tx, "sampling_frequency", &actual_fs),
                      "Read sample rate");
      check_iio_error(
          iio_channel_attr_read_longlong(phy_tx, "rf_bandwidth", &actual_bw),
          "Read RF bandwidth");
      if (actual_fs != sample_rate)
        throw std::runtime_error("Effective sample rate differs from 2400000: " +
                                 std::to_string(actual_fs));
      std::cout << "Effective sample rate: " << actual_fs
                << " Hz; RF filter: " << actual_bw << " Hz\n";

      for (unsigned n = 0; n < iio_device_get_channels_count(tx); ++n) {
        auto *channel = iio_device_get_channel(tx, n);
        if (iio_channel_is_scan_element(channel))
          iio_channel_disable(channel);
      }
      tx_i_ = iio_device_find_channel(tx, "voltage0", true);
      tx_q_ = iio_device_find_channel(tx, "voltage1", true);
      if (!tx_i_ || !tx_q_)
        throw std::runtime_error("TX IQ channels not found");
      scale_i_ = channel_scale(tx_i_);
      scale_q_ = channel_scale(tx_q_);
      iio_channel_enable(tx_i_);
      iio_channel_enable(tx_q_);
      buffer_ = iio_device_create_buffer(tx, buffer_samples, false);
      if (!buffer_)
        throw std::runtime_error("Failed to create TX buffer");
    } catch (...) {
      cleanup();
      throw;
    }
  }
  ~PlutoTransmitter() { cleanup(); }
  PlutoTransmitter(const PlutoTransmitter &) = delete;
  PlutoTransmitter &operator=(const PlutoTransmitter &) = delete;

  void transmit_buffer(const std::vector<IQ> &samples) {
    if (samples.size() != buffer_samples)
      throw std::runtime_error("The buffer must be filled exactly");
    auto *begin = static_cast<char *>(iio_buffer_start(buffer_));
    auto *end = static_cast<char *>(iio_buffer_end(buffer_));
    auto *ip = static_cast<char *>(iio_buffer_first(buffer_, tx_i_));
    auto *qp = static_cast<char *>(iio_buffer_first(buffer_, tx_q_));
    const auto step = iio_buffer_step(buffer_);
    if (!begin || !end || !ip || !qp || step < 4 || end < begin ||
        static_cast<size_t>(end - begin) / static_cast<size_t>(step) !=
            samples.size() ||
        ip < begin || qp < begin || ip >= end || qp >= end ||
        (samples.size() - 1) * static_cast<size_t>(step) + sizeof(int16_t) >
            static_cast<size_t>(end - ip) ||
        (samples.size() - 1) * static_cast<size_t>(step) + sizeof(int16_t) >
            static_cast<size_t>(end - qp))
      throw std::runtime_error("Invalid TX buffer format");
    for (size_t n = 0; n < samples.size(); ++n) {
      const int16_t i = quantize(samples[n].i, scale_i_);
      const int16_t q = quantize(samples[n].q, scale_q_);
      iio_channel_convert_inverse(tx_i_, ip + n * step, &i);
      iio_channel_convert_inverse(tx_q_, qp + n * step, &q);
    }
    const auto pushed = iio_buffer_push(buffer_);
    check_iio_error(pushed, "TX transfer");
    if (pushed != end - begin)
      throw std::runtime_error("Incomplete TX transfer");
  }

private:
  static int channel_scale(const iio_channel *channel) {
    const auto *format = iio_channel_get_data_format(channel);
    if (!format || format->length != 16 || !format->is_signed ||
        format->bits < 2 || format->bits > 16 || format->repeat > 1 ||
        format->bits + format->shift > format->length)
      throw std::runtime_error("Expected signed 16-bit IQ format");
    return (1 << (format->bits - 1)) - 1;
  }
  static int16_t quantize(float value, int scale) {
    return static_cast<int16_t>(
        std::lrint(std::clamp(value, -1.0f, 1.0f) * scale));
  }
  void cleanup() noexcept {
    if (buffer_) {
      iio_buffer_destroy(buffer_);
      buffer_ = nullptr;
    }
    if (tx_i_) {
      iio_channel_disable(tx_i_);
      tx_i_ = nullptr;
    }
    if (tx_q_) {
      iio_channel_disable(tx_q_);
      tx_q_ = nullptr;
    }
    if (ctx_) {
      iio_context_destroy(ctx_);
      ctx_ = nullptr;
    }
  }
  iio_context *ctx_ = nullptr;
  iio_channel *tx_i_ = nullptr;
  iio_channel *tx_q_ = nullptr;
  iio_buffer *buffer_ = nullptr;
  int scale_i_ = 0, scale_q_ = 0;
};
#endif

uint64_t parse_integer(const std::string &text, uint64_t maximum) {
  if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
    throw std::runtime_error("Expected an unsigned integer: " + text);
  size_t used = 0;
  const auto value = std::stoull(text, &used);
  if (used != text.size() || value > maximum)
    throw std::runtime_error("Integer out of range: " + text);
  return value;
}
}

#ifndef OQPSK_NO_MAIN
int main(int argc, char **argv) {
  using namespace oqpsk;
  std::signal(SIGINT, signal_handler);
  std::signal(SIGTERM, signal_handler);
  try {
    const std::string mode = argc > 1 ? argv[1] : "--help";
    if (mode == "--help") {
      std::cout << "Usage:\n  " << argv[0] << " --tx [URI]\n  " << argv[0]
                << " --dump file.cf32 sample_count [seed]\n"
                << "Example: --dump oqpsk.cf32 2400000 12345\n"
                << "Interleaved float32 IQ, native byte order; Fs = 2400000 Hz.\n";
      return 0;
    }
    if ((mode == "--tx" && argc > 3) ||
        (mode == "--dump" && argc != 4 && argc != 5) ||
        (mode != "--tx" && mode != "--dump"))
      throw std::runtime_error("Invalid arguments; use --help");
#ifdef OQPSK_OFFLINE
    if (mode == "--tx")
      throw std::runtime_error("Recompile with libiio to transmit");
#endif
    const uint32_t seed =
        argc == 5 ? static_cast<uint32_t>(parse_integer(
                        argv[4], std::numeric_limits<uint32_t>::max()))
                  : std::random_device{}();
    Modulator modulator(seed);
    modulator.prime();
    std::vector<IQ> samples(buffer_samples);
    std::cout << "OQPSK : 10000 sym/s, 20000 bit/s, RRC beta=0.25, span=16\n"
              << "Theoretical bandwidth: 12500 Hz (+/-6250 Hz), sample rate: 2400000 Hz\n"
              << "Payload: concatenated blocks of 1024 random bits; seed="
              << seed << "\n";
    if (mode == "--dump") {
      const uint64_t count = parse_integer(
          argv[3], std::numeric_limits<uint64_t>::max() / sizeof(IQ));
      if (!count)
        throw std::runtime_error("Sample count must be nonzero");
      std::ofstream out(argv[2], std::ios::binary | std::ios::trunc);
      if (!out)
        throw std::runtime_error("Failed to open IQ output file");
      uint64_t written = 0;
      while (written < count && !stop_requested) {
        const auto n = static_cast<size_t>(
            std::min<uint64_t>(samples.size(), count - written));
        modulator.fill(samples.data(), n);
        out.write(reinterpret_cast<const char *>(samples.data()),
                  static_cast<std::streamsize>(n * sizeof(IQ)));
        if (!out)
          throw std::runtime_error("Failed to write IQ output");
        written += n;
      }
      out.close();
      if (!out)
        throw std::runtime_error("Failed to close IQ output file");
      std::cout << written << " IQ samples written\n";
    }
#ifndef OQPSK_OFFLINE
    else {
      PlutoTransmitter transmitter(argc == 3 ? argv[2] : "usb:");
      std::cout << "TX : " << tx_frequency << " Hz, gain : " << tx_gain_db
                << " dB. Press Ctrl-C to stop.\n";
      while (!stop_requested) {
        modulator.fill(samples.data(), samples.size());
        if (!stop_requested)
          transmitter.transmit_buffer(samples);
      }
    }
#endif
  } catch (const std::exception &error) {
    if (stop_requested) {
      std::cerr << "Stop requested.\n";
      return 0;
    }
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
#endif
