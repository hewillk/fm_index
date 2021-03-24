## Biomodern.FmIndex
`biomodern::FMIndex` is specialization of FM-Index for genome sequence which implemented using C++20. The core data structure of `FMIndex` are following:

```cpp
biomodern::DibitVector bwt; // a compression vector which store the bwt.
vector<array<uint32_t, 4>> occ1;
vector<array<uint8_t, 4>> occ2; // a hierarchical sampled occurrence table.
vector<uint32_t> sa; // a sampled suffix array.
vector<uint32_t> lookup; // a lookup table for fixed suffix query.
```

The sampling parameters are those three:
- `occ_intv`: occ sampling interval, default value is `16`.
- `sa_intv`: sa sampling interval, default value is `1`.
- `lookup_len`: the length of fixed suffix for lookup, default value is `14`.

The memory usage for `FMIndex` in run time is affect by above parameters. Take a `3.1Gb` human genome for example, the memory occupation can be calculate as follwing:
- bwt string: `3.1Gb / 4 = 0.775 Gb`.
- hierarchical occurrence table: L1 occ occupy fixed `3.1Gb * 16 / 256 = 0.194Gb` plus L2 occ occupy `3.1Gb * 4 / occ_intv(16) = 0.775 Gb`.
- suffix array: `3.1Gb * 4 / sa_intv(1) = 12Gb`. Noticed that in default mode we dont sampling suffix value since this can reduce frequently memory allocation and intense computation when occurs massive query.
- lookup table: `4^lookup_len(14) * 4 / 1024^3 = 1Gb`.

So the default total memory occupation for `3.1Gb` human genome is `0.775 Gb + 0.194Gb + 0.775Gb + 12Gb + 1Gb = 14.744Gb`.

## Compilers
- GCC 10.2

## Usage
The FMIndex class provide very simple interface for ultra fast exact match:
```cpp
using istring_view = basic_string_view<int8_t>;

// Utility for compute the suffix array from istring
static vector<uint32_t> FMIndex::get_sa(istring_view ref);

// Ctor. initialize the core data structure from reference istring. 
// Noticed that this Ctor will call FMIndex::get_sa(...) to generate suffix array.
FMIndex(istring_view ref);

// Get begin and end index of the uncompressed suffix array for the input seed.
// Difference between begin and end is the occurrence count in reference.
// The third value is the remaining prefix length for the seed. 
// Notice that when the stop_cnt is -1, this value always be 0.
// The stop_cnt can be using to early stop when occurrence count is not greater than the value.
// Set to -1 to forbid early stop.
tuple<uint32_t, uint32_t, uint32_t> get_range(istring_view seed, uint32_t stop_cnt);

// Compute the suffix array value according to the begin and end. 
// If sa_intv is 1, this can be done at O(1).
vector<uint32_t> get_offsets(uint32_t, uint32_t);

// Utility for serialization. 
// The binary binary archive file size is same as the memory occupation in run time.
void save(ofstream&);
void load(ifstream&);
```

## Example
```cpp
#include <iostream>
#include <experimental/random>
#include "fm_index.hpp"
  
int main() {
  using namespace biomodern::utility;
  using namespace std::string_literals;

  {
    auto ref = istring{};
    ref.reserve(3137454505);
    auto fin = std::ifstream{"hs37d5.fa"};
    for (auto line = ""s; std::getline(fin, line);) {
      if (line.front() == '>') {
        std::cout << line << "\n";
        continue;
      }
      auto subref = Codec::to_istring(line);
      // fm-index only support four characters so we need change 'N' to 'ACGT'
      for (auto& c : subref) if (c == 4) c = std::experimental::randint(0, 3);
      ref += subref;
    }
    const auto fmi = biomodern::FMIndex{ref};
    auto fout = std::ofstream{"hs37d5.fmi", std::ios::binary};
    fmi.save(fout);
  }
  
  auto fmi = biomodern::FMIndex{};
  {
    auto fin = std::ifstream{"hs37d5.fmi", std::ios::binary};
    fmi.load(fin);
  }
  
  for (auto seed = istring{}; std::cin >> seed;) {
    const auto [beg, end, offset] = fmi.get_range(seed, 0);
    std::cout << "seed: " << seed << "\n";
    std::cout << "seed offset: " << offset << "\n";
    std::cout << "occurrence: " << end - beg << "\n";
    const auto offsets = fmi.get_offsets(beg, end);
    std::cout << "ref offsets: ";
    for (const auto offset : offsets) std::cout << offset << " ";
      std::cout << "\n";
  }
}
```

