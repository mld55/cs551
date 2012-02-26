#define DEBUG_TYPE "cs551a1"
#include "llvm/Module.h"
#include "llvm/Pass.h"
#include "llvm/PassSupport.h"
#include "llvm/Type.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/Support/Debug.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/ADT/Statistic.h"

using namespace llvm;

STATISTIC(CS551A1Count, "Counts number of loops");
STATISTIC(CS551A1Depth, "Counts depth of loops");

namespace {
    struct CS551A1 : public LoopPass
    {
        static char ID;
        CS551A1() : LoopPass(ID) {}

        virtual bool runOnLoop(Loop *LP, LPPassManager &LPM);
        virtual void print(raw_ostream &O, const Module *M) const;
        private:
        LoopInfo *LI;
        Loop *LP;
    };
}

// I do not understand what this does;
// it is just copy-and-pasted from the example :-(
char CS551A1::ID = 0;
static RegisterPass<CS551A1> X("cs551a1", "Loop Count and Depth");

bool CS551A1::runOnLoop(Loop *LP, LPPassManager &LPM)
{
    this->LP = LP;
    // return if we modified the graph
    bool ModifiedP = false;
    ++CS551A1Count;
    CS551A1Depth += LP->getLoopDepth();
    return ModifiedP;
}

void CS551A1::print(raw_ostream &O, const Module *M) const
{
    O << "hello, I am a CS551a1 module\n";
}
