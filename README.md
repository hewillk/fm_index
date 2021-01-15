## biomodern-fm-index
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
- `lookup_len`: the length of fixed suffix for lookup, default value is `13`.

The memory usage for `FMIndex` in run time is affect by above parameters. Take a `3.1Gb` human genome for example, the memory occupation can be calculate as follwing:
- bwt string: `3.1Gb / 4 = 0.775 Gb`.
- hierarchical occurrence table: L1 occ occupy fixed `3.1Gb * 16 / 256 = 0.194Gb` plus L2 occ occupy `3.1Gb * 4 / occ_intv(16) = 0.775 Gb`.
- suffix array: `3.1Gb * 4 / sa_intv(1) = 12Gb`. Noticed that in default mode we dont sampling suffix value since this can reduce frequently memory allocation and intense computation when occurs massive query.
- lookup table: `4^lookup_len(13) * 4 / 1024^3 = 0.25Gb`.

So the default total memory occupation of `FMIndex` for `3.1Gb` human genome is `0.775 Gb + 0.194Gb + 0.775Gb + 12Gb + 0.25Gb = 13.994Gb`.

## Usage
The FMIndex class provide very simple interface for ultra fast exact match:
```cpp
using istring_view = basic_string_view<int8_t>;

// Utility for compute the suffix array from istring
static vector<uint32_t> FMIndex::get_sa(istring_view ref);

// Ctor. initialize the core data structure from reference istring. 
// Noticed that this Ctor will call FMIndex::get_sa(...) to generate suffix array.
FMIndex(istring_view ref);

// Return begin and end index of the uncompressed suffix array for the input query.
// Difference between begin and end is the occurrence count in reference.
// The stop_cnt can be using to early stop when occurrence count is not greater than the value.
pair<uint32_t, uint32_t> get_range(istring_view seed, uint32_t stop_cnt);

// Compute the suffix array value according to the begin and end. 
// If sa_intv is 1, this can be done at O(1).
vector<uint32_t> get_offsets(uint32_t, uint32_t);

// Utility for serialization. 
// The binary binary archive file size is same as the memory occupation in run time.
void save(ofstream&);
void load(ifstream&);
```

