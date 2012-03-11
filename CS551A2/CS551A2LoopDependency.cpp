#define DEBUG_TYPE "cs551a1"

#include "llvm/Module.h"
#include "llvm/Pass.h"
#include "llvm/PassSupport.h"
#include "llvm/Type.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/Analysis/LoopDependenceAnalysis.h"
#include "llvm/Support/Debug.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/ADT/Statistic.h"

using namespace llvm;

namespace {
    struct CS551A2 : public LoopDependenceAnalysis
    {
        static char ID;
        CS551A2() : LoopDependenceAnalysis() {}
        virtual bool runOnLoop(Loop *LP, LPPassManager &LPM);
        virtual void print(raw_ostream &O, const Module *M) const;
        private:
        LoopInfo *LI;
        Loop *LP;
    };
}

// I do not understand what this does;
// it is just copy-and-pasted from the example :-(
char CS551A2::ID = 0;
static RegisterPass<CS551A2> X("cs551a2", "Loop Dependency Magick");

bool CS551A2::runOnLoop(Loop *LP, LPPassManager &LPM)
{
    this->LP = LP;
    errs() << "Hello from MyLDA\n";
    // return if we modified the graph
    bool ModifiedP = false;
    return ModifiedP;
}

void CS551A2::print(raw_ostream &O, const Module *M) const
{
    O << "hello, I am a CS551a2 module\n";
}
