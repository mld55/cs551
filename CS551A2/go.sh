#! /bin/sh
PATH=$HOME/Downloads/xcode-llvm/bin:$HOME/Downloads/xcode-llvm/bin/Debug:$PATH
ll_libdir=$HOME/Downloads/xcode-llvm/lib/Debug
ll_libname=LLVMCS551A2.dylib 
# =None       -   disable debug output
# =Arguments  -   print pass arguments to pass to 'opt'
# =Structure  -   print pass structure before run()
# =Executions -   print pass name before it is executed
# =Details    -   print pass details when it is executed
#    -lda \
ulimit -c unlimited
# basicaa comes from ``https://groups.google.com/forum/?fromgroups&hl=en#!topic/llvm-dev/rpO2_iNjk74''
opt_argv="$opt_argv -basicaa"
# opt_argv="$opt_argv -licm"
# opt_argv="$opt_argv -aa-eval"
opt_argv="$opt_argv -print-alias-sets"
# opt_argv="$opt_argv -indvars -enable-iv-rewrite"
case "$1" in
    mine)
    opt_argv="$opt_argv -load ${ll_libdir}/${ll_libname} -cs551a2"
    ;;
    lda)
    opt_argv="$opt_argv -lda"
    ;;
esac
# opt_argv="$opt_argv --debug-pass=Details"
opt_argv="$opt_argv -debug"
cat $2 | opt ${opt_argv} -disable-output
