#ifndef BIOMODERN__FM_INDEX_HPP_
#define BIOMODERN__FM_INDEX_HPP_

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <execution>
#include <iostream>
#include <numeric>
#include <span>
#include <string>
#include <tuple>
#include <vector>

#include "dibit_vector/dibit_vector.hpp"
#include "istring/istring.hpp"
#include "serializer/serializer.hpp"

namespace biomodern {

using namespace std::chrono;
using namespace utility;

class FMIndex {
 public:
  using size_type = std::uint32_t;
  using char_type = std::int8_t;

  constexpr static auto occ_intv = 16u;
  constexpr static auto sa_intv = 1u;
  constexpr static auto lookup_len = 13u;

  constexpr static auto sort_len = 256u;
  constexpr static auto occ1_intv = 256u;
  constexpr static auto occ2_intv = occ_intv;

  static auto get_sa(istring_view ref) {
    auto sa = std::vector<size_type>(ref.size() + 1);
    std::iota(sa.begin(), sa.end(), 0);
    std::cout << "sa sort start...\n";
    const auto start = high_resolution_clock::now();
    std::sort(std::execution::par_unseq, sa.begin(), sa.end(),
              [ref](const auto i, const auto j) {
                return ref.substr(i, sort_len) < ref.substr(j, sort_len);
              });
    const auto end = high_resolution_clock::now();
    const auto dur = duration_cast<seconds>(end - start);
    std::cout << "sorting time: " << dur.count() << " s.\n";
    return sa;
  }

 private:
  // core data structure
  DibitVector<std::uint8_t> bwt_;
  std::pair<std::vector<std::array<size_type, 4>>,
            std::vector<std::array<std::uint8_t, 4>>>
      occ_;
  std::vector<size_type> sa_;
  std::array<size_type, 4> cnt_{};
  size_type pri_{};
  std::vector<size_type> lookup_;

  constexpr static auto cnt_table = [] {
    std::array<std::array<std::uint8_t, 4>, 256> cnt_table{};
    for (auto i = size_type{}; i < cnt_table.size(); i++)
      for (auto shift = size_type{}; shift < 8; shift += 2)
        cnt_table[i][i >> shift & 3u]++;
    return cnt_table;
  }();

  auto compute_occ(char_type c, size_type i) const {
    const auto occ1_beg = i / occ1_intv;
    const auto occ2_beg = i / occ2_intv;
    auto beg = occ2_beg * occ2_intv;
    auto cnt = size_type{};
    const auto pass_pri = c == 0 && beg <= pri_ && pri_ < i;
    const auto run = (i - beg) / 4;
    for (auto j = size_type{}; j < run; j++) {
      cnt += cnt_table[bwt_.data()[beg / 4]][c];
      beg += 4;
    }
    for (; beg < i; beg++)
      if (bwt_[beg] == c) cnt++;
    return occ_.first[occ1_beg][c] + occ_.second[occ2_beg][c] + cnt - pass_pri;
  }

  auto lf(char_type c, size_type i) const {
    return cnt_[c] + compute_occ(c, i);
  };

  auto compute_sa(size_type i) const {
    auto cnt = size_type{};
    while (i % sa_intv && i != pri_) {
      i = lf(bwt_[i], i);
      cnt++;
    }
    return i != pri_ ? sa_[i / sa_intv] + cnt : cnt;
  }

  auto compute_range(istring_view seed, size_type beg, size_type end,
                     size_type stop_upper) const {
    while (!seed.empty()) {
      if (end - beg < stop_upper) break;
      beg = lf(seed.back(), beg);
      end = lf(seed.back(), end);
      seed.remove_suffix(1);
    }
    return std::tuple{beg, end, seed.size()};
  }

  auto compute_lockup() {
    constexpr auto lookup_size = size_type{1} << lookup_len * 2;
    lookup_.reserve(lookup_size + 1);
    std::cout << "computing " << lookup_size << " suffix for for lockup...\n";
    const auto start = high_resolution_clock::now();
    std::cout << "compute start...\n";
    for (auto i = size_type{}; i < lookup_size; i++) {
      const auto key = Codec::rhash(i, lookup_len);
      const auto [beg, end, offset] = compute_range(key, 0, bwt_.size(), 0);
      lookup_.push_back(beg);
    }
    lookup_.push_back(bwt_.size());
    const auto end = high_resolution_clock::now();
    const auto dur = duration_cast<seconds>(end - start);
    std::cout << "computing time: " << dur.count() << " s.\n";
  }

 public:
  FMIndex() = default;
  FMIndex(istring_view ref) {
    std::cout << "validate ref...\n";
    assert(std::ranges::all_of(ref, [](auto c) { return c >= 0 && c <= 3; }));

    auto ori_sa = get_sa(ref);
    bwt_.reserve(ori_sa.size());
    auto& [occ1, occ2] = occ_;
    occ1.reserve(ori_sa.size() / occ1_intv + 1);
    occ2.reserve(ori_sa.size() / occ2_intv + 1);
    sa_.reserve(ori_sa.size() / sa_intv + 1);
    auto& cnt1 = cnt_;
    auto cnt2 = std::array<std::uint8_t, 4>{};
    occ1.push_back(cnt1);
    occ2.push_back(cnt2);
    for (auto i = size_type{}; i < ori_sa.size(); i++) {
      const auto sa_v = ori_sa[i];
      if (sa_v != 0) {
        const auto c = ref[sa_v - 1];
        bwt_.push_back(c);
        cnt1[c]++;
        cnt2[c]++;
      } else {
        bwt_.push_back(0);
        pri_ = i;
      }
      if ((i + 1) % occ2_intv == 0) {
        if ((i + 1) % occ1_intv == 0) {
          occ1.push_back(cnt1);
          cnt2 = {};
        }
        occ2.push_back(cnt2);
      }
      if constexpr (sa_intv != 1) {
        if (i % sa_intv == 0) sa_.push_back(sa_v);
      }
    }
    auto sum = size_type{1};
    for (auto& x : cnt1) {
      sum += x;
      x = sum - x;
    }
    if constexpr (sa_intv == 1) sa_.swap(ori_sa);
    compute_lockup();
  }

  auto get_offsets(size_type beg, size_type end) const {
    if constexpr (sa_intv == 1)
      return std::span{&sa_[beg], end - beg};
    else {
      auto offsets = std::vector<size_type>{};
      offsets.reserve(end - beg);
      for (auto i = beg; i < end; i++) offsets.push_back(compute_sa(i));
      std::ranges::sort(offsets);
      return offsets;
    }
  }

  auto get_range(istring_view seed, size_type beg, size_type end,
                 size_type stop_cnt) const {
    if (end == beg || seed.empty()) return std::tuple{beg, end, 0ul};
    return compute_range(seed, beg, end, stop_cnt + 1);
  }

  auto get_range(istring_view seed, size_type stop_cnt) const {
    auto beg = size_type{};
    auto end = bwt_.size();
    if (seed.size() >= lookup_len) {
      const auto key = Codec::hash(seed.substr(seed.size() - lookup_len));
      beg = lookup_[key];
      end = lookup_[key + 1];
      seed.remove_suffix(lookup_len);
    }
    return get_range(seed, beg, end, stop_cnt);
  }

  auto save(std::ofstream& fout) const {
    fout.write(reinterpret_cast<const char*>(&cnt_), sizeof(cnt_));
    fout.write(reinterpret_cast<const char*>(&pri_), sizeof(pri_));
    std::cout << "save bwt...\n";
    Serializer::save(fout, bwt_);
    std::cout << "save occ...\n";
    Serializer::save(fout, occ_.first);
    Serializer::save(fout, occ_.second);
    std::cout << "save sa...\n";
    Serializer::save(fout, sa_);
    std::cout << "save lockup...\n";
    Serializer::save(fout, lookup_);
  }

  auto load(std::ifstream& fin) {
    fin.read(reinterpret_cast<char*>(&cnt_), sizeof(cnt_));
    fin.read(reinterpret_cast<char*>(&pri_), sizeof(pri_));
    std::cout << "load bwt...\n";
    Serializer::load(fin, bwt_);
    std::cout << "load occ...\n";
    Serializer::load(fin, occ_.first);
    Serializer::load(fin, occ_.second);
    std::cout << "load sa...\n";
    Serializer::load(fin, sa_);
    std::cout << "load lookup...\n";
    Serializer::load(fin, lookup_);
    assert(fin.peek() == EOF);
  }

  bool operator==(const FMIndex& other) const = default;
};

}  // namespace biomodern

#endif
