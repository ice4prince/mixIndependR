## Test environments

* x86_64-w64-mingw32 (64-bit), R 4.0.3
* x86_64-apple-darwin17.0 (64-bit), R 4.0.3


* Rhub: Platform: Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Rhub: Platform: Fedora Linux, R-devel, clang, gfortran


## R CMD check results


0 errors ✓ | 0 warnings ✓ | 0 notes ✓

## rhub check results  

  Maintainer: 'Bing Song <bingsong@my.unthsc.edu>'
  
  Days since last update: 3

0 errors ✓ | 0 warnings ✓ | 1 note x

I appologize for the constant update. I was in a rush to submit revised manuscript to journal and failed to find the bugs on documentation of "read_vcf_gt" and "splitGenotypes". I forgot to export them to the NAMESPACE. Since they are not used in vignettes or examples, these bugs had been overlooked. But in real application, it matters. This update corrected this.

