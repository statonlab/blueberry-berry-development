## errors when running LoReAn
the program break at the trinity part. 
```text
salmon: error while loading shared libraries: liblzma.so.0: cannot open shared object file: No such file or directory
```
Search about the error, most answers said 'liblzma.so.0 should be supplied and owned by the xz-libs package.'. So get out of conda environment and try again.
