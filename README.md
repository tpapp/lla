# Lisp Linear Algebra — a linear algebra library for Common Lisp

[![Project Status: Unsupported – The project has reached a stable, usable state but the author(s) have ceased all work on it. A new maintainer may be desired.](http://www.repostatus.org/badges/latest/unsupported.svg)](http://www.repostatus.org/#unsupported) This library is [**unsupported**](https://tpapp.github.io/post/orphaned-lisp-libraries/).

LLA is a high-level Common Lisp library built on on [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/), but providing a much more abstract interface with the purpose of freeing the user from low-level concerns and reducing the number of bugs in numerical code.

Documentation is mostly in docstrings at the moment, but I plan to write a decent tutorial at some point. In the meantime, please look at the unit tests.

## Objectives

- High-level, user friendly interface that hides the details.

    `(solve a b)` should return $X$, from $AX=B$, regardless of whether $A$ is a dense matrix, an $LU$ decomposition, or something else; similarly, $X$ should be a vector/matrix when $B$ is. Users should not need to memorize names like `DGESV`, especially when CLOS makes it so easy to deal with these things. Also, you don't need to make sure that all arguments are of the same type (eg complex-double): LLA will find the common supertype for elements and convert if necessary.

- Stay in Lisp and expose the innards of the objects as much as possible.

    LLA aims to take advantage of CL's high level facilities such as CLOS and memory management. Data is kept in Lisp arrays instead of foreign arrays, so you can access it directly using `aref` etc. You also benefit from garbage collection and all the clever stuff that comes with the GC. If you need more memory, just increase the heap size.

- Keeping it simple.

    Currently, LLA sources amount to less than 3000 lines of code (not including tests). The small size should make maintainance easier and bugs more rare (hopefully).

- Speed is important, but reliability comes first.

    Only optimize when necessary, and do extensive testing afterwards. Most of the speed comes from your LAPACK library anyway --- most linear algebra operations are $O(N^\alpha)$ with $\alpha > 1$, frequently $\alpha > 2$. That said, copying to memory is optimized, and in the long run LLA should make use of your implementation's array pinning features (when available). *Currently, direct array sharing is disabled, it will be re-enabeld in the near future*.

# Configuration

Certain features of LLA can be configured before loading using the plist `*lla-configuration*` in the `CL-USER` package (for example, on SBCL you would do it in your `~/.sbclrc`). The following properties are supported:

- `:libraries`

    A list of objects, passed directly to `cffi:load-foreign-library`. You can use strings, paths, or even symbols if you have defined these libraries using `cffi:define-foreign-library`. If you don't define this, a reasonable platform-dependent default will be used. See the next section for details.

- `:int64`

    This makes LLA use 64-bit integers for array dimensions, pivot indices and other integer values passed to BLAS/LAPACK functions. **Only use this if you are sure that your libraries have been compiled with 64-bit integers**. The fact that you have a 64-bit platform does not necessarily mean that this is the case, in fact, it is still quite rare. Unless told otherwise, LLA expectes BLAS/LAPACK to use the [(L)LP64 model](http://en.wikipedia.org/wiki/64-bit#64-bit_data_models) for integers -- that is to say, integer types in Fortran are 32 bit.

- `:efficiency-warnings`

    Enable the **possibility** of efficiency warnings at compile time. You still have to set the appropriate flags, but without this option, they won't even be checked. There are two properties you can set: `:array-type` and `:array-conversion`. The first warns whenever an array has to be walked elementwise to determine its type, the second when some arrays need to be converted to a common type.

    Example:

    ```lisp
    (defparameter cl-user:*lla-configuration*
      '(:efficiency-warnings (:array-type :array-conversion)))
    ```

    before loading LLA, and

    ```lisp
    (let ((lla:*lla-efficiency-warning-array-type* t)
          (lla:*lla-efficiency-warning-array-conversion* t))
       (code that you want to check))
    ```

## Libraries

## Dependencies and configuration

LLA needs BLAS and LAPACK shared libraries to work. When it comes to
loading libraries, LLA tries to pick a sensible default for each
platform, but in case it fails, you need to tell LLA where the libraries
are before loading.

You can do this by putting something like this in your startup script (eg `~/.sbclrc`, the symbol needs to be in the package `cl-user`):

```lisp
(defvar *lla-configuration*
  '(:libraries ("/usr/lib/atlas-base/atlas/libblas.so.3gf"
                "/usr/lib/atlas-base/libatlas.so.3gf")))
```

## Debian

On Debian-based distributions, it is very likely that LLA will work out of the box if you just install ATLAS, eg

```sh
apt-get install libatlas3gf-base
```

However, you may want to build a version optimized for your architecture.

### Building ATLAS on Debian

Prepare the build (as root):

```sh
apt-get build-dep atlas
apt-get install fakeroot devscripts
cpufreq-set -g performance -c 0   # do this for all CPUs
```

Then as a regular user,

```sh
apt-get source atlas
cd atlas-[fill in your version here]/
fakeroot debian/rules custom
```

Then install the .deb files that were created.

### Selecting the right linear algebra library

```sh
update-alternatives --config libblas.so.3
update-alternatives --config liblapack.so.3
```

### Intel MKL on Linux

In `/etc/ld.so.conf.d/`, create a file that contains the paths, eg

```
/opt/intel/mkl/lib/intel64
/opt/intel/composerxe/lib/intel64
```

Then the configuration

```lisp
(defvar *lla-configuration*
  '("libgomp.so.1" "libiomp5.so" "libmkl_rt" "libpthread.so.0" "libpthread"))
```

should work.

## Acknowledgements

LLA was inspired by packages written by AJ Rossini, Rif, Mark Hoemmen and others. I have borrowed code (whenever allowed by their licenses) and ideas freely from all of them.

Gábor Melis made substantial contributions to the library, especially the low-level pinning interface and the destructive BLAS routines.

## Suggested editor settings for code contributions

No line breaks in (doc)strings, otherwise try to keep it within 80 columns. Remove trailing whitespace. 'modern' coding style. Suggested Emacs snippet:

```emacs-lisp
(set-fill-column 9999)
(font-lock-add-keywords nil
                        '(("\\<\\(FIXME\\|TODO\\|QUESTION\\|NOTE\\)"
                        1 font-lock-warning-face t)))
(setq show-trailing-whitespace t)
(add-hook 'write-file-hooks
          '(lambda()
             (save-excursion
               (delete-trailing-whitespace))
             nil))
(visual-line-mode 1)
(setq slime-net-coding-system 'utf-8-unix)
(setq lisp-lambda-list-keyword-parameter-alignment t)
(setq lisp-lambda-list-keyword-alignment t)
(setq common-lisp-style-default 'modern)
```

## Things to do (roughly in order of priority)

-   write optimized pinning interfaces, especially ECL
-   write documentation (probably w/ [docudown](http://common-lisp.net/project/docudown/), decide)
-   write more tests (especially randomized ones, develop macros for
    that)
-   write a tutorial
