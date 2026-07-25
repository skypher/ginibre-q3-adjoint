# Reconstruct the exact `O_13` low-dimensional replay archive

```sh
cat su2_o13_lowdim_exact_replay.tar.gz.b64.part* \
  | base64 -d > su2_o13_lowdim_exact_replay.tar.gz
sha256sum -c SHA256SUMS

tar -xzf su2_o13_lowdim_exact_replay.tar.gz
./replay_o13_lowdim/replay.sh
```
