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
- hierarchical occurrence table: L1 occ occupy fixed `3.1Gb * 16 / 256 = 0.194Gb` plus variable L2 occ occupy `3.1Gb * 4 / occ_intv(16) = 0.775 Gb`.
- suffix array: `3.1Gb * 4 / sa_intv(1) = 12Gb`. Noticed that in default mode we dont sampling suffix value which can reduce frequently memory allocation and intense computation when occurs massive query.
- lookup table: `4^lookup_len(13) * 4 / 1024^3 = 0.25Gb`.

So the default total memory occupation of `FMIndex` for `3.1Gb` human genome is `0.775 Gb + 0.194Gb + 0.775Gb + 12Gb + 0.25Gb = 13.994Gb`.




