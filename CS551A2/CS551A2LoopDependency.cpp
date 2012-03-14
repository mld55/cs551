#define DEBUG_TYPE "cs551a1"

#include "llvm/Module.h"
#include "llvm/Pass.h"
#include "llvm/PassSupport.h"
#include "llvm/Type.h"

#include "llvm/ADT/DenseSet.h"
#include "llvm/ADT/FoldingSet.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/Statistic.h"

#include "llvm/Analysis/AliasAnalysis.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/Analysis/ScalarEvolution.h"

#include "llvm/Support/Allocator.h"
#include "llvm/Support/Debug.h"
#include "llvm/Support/raw_ostream.h"
// IMPL {{{
#include "llvm/Instructions.h"
#include "llvm/Operator.h"
#include "llvm/Analysis/ScalarEvolutionExpressions.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/Assembly/Writer.h"
#include "llvm/Support/ErrorHandling.h"
#include "llvm/Target/TargetData.h"
// IMPL }}}

using namespace llvm;

namespace
{
/*
    class Assignment2
    {
        enum DirectionType {    
            LT, EQ, GT, STAR
        };
        struct DirectionV
        {
            std::vector<DirectionType> dir;
        };
        struct DistanceV
        {
            std::vector<int> distance;
        };
        bool isSeparable(Subscript* S);
        bool isZIV(Subscript* S);
        bool isSIV(Subscript* S);
        bool isMIV(Subscript* S);
        DirectionV *getDirectionVector(Subscript* S);
        DistanceV *getDistanceVector(Subscript* S);
        void partition(std::vector<Subscript*> S,
                       std::vector< std::vector<Subscript*> > *P,
                       int *nP);
        void delta_test(std::vector<Subscript*> &subscripts,
                        std::vector<DirectionV*> *DVset,
                        std::vector<DistanceV*> *dVset);
    };
*/
    class CS551A2 : public LoopPass
    {
    public:
        struct Affine
        {
            /// contains the loop index coefficient, or 1
            int coefficient;
            const Value *index;
            /// contains the constant added to the operand above, or 0
            int constant;
        };
        struct Subscript
        {
            const Affine *A;
            const Affine *B;
        };

    public:
        // BEGIN COPY-PASTE {{{
        /// TODO: doc
        enum DependenceResult { Independent = 0, Dependent = 1, Unknown = 2 };

        /// DependencePair - Represents a data dependence relation between to memory
        /// reference instructions.
        struct DependencePair : public FastFoldingSetNode
        {
            Value *A;
            Value *B;
            DependenceResult Result;
            DependencePair(const FoldingSetNodeID &ID, Value *a, Value *b)
                : FastFoldingSetNode(ID), A(a), B(b), Result(Unknown)
                {}
        };

        /// findOrInsertDependencePair - Return true if a DependencePair for the
        /// given Values already exists, false if a new DependencePair had to be
        /// created. The third argument is set to the pair found or created.
        bool findOrInsertDependencePair(Value*, Value*, DependencePair*&);

        /// getLoops - Collect all loops of the loop nest L in which
        /// a given SCEV is variant.
        void getLoops(const SCEV*, DenseSet<const Loop*>*) const;

        /// isLoopInvariant - True if a given SCEV is invariant in all loops of the
        /// loop nest starting at the innermost loop L.
        bool isLoopInvariant(const SCEV*) const;

        /// isAffine - An SCEV is affine with respect to the loop nest starting at
        /// the innermost loop L if it is of the form A+B*X where A, B are invariant
        /// in the loop nest and X is a induction variable in the loop nest.
        bool isAffine(const SCEV*) const;

        Affine* getAffine(const Value *);

        /// TODO: doc
        bool isZIVPair(const SCEV*, const SCEV*) const;
        bool isSIVPair(const SCEV*, const SCEV*) const;
        DependenceResult analyseZIV(const SCEV*, const SCEV*) const;
        DependenceResult analyseSIV(const SCEV*, const SCEV*) const;
        DependenceResult analyseMIV(const SCEV*, const SCEV*) const;
        DependenceResult analyseSubscript(const SCEV*, const SCEV*) const;
        DependenceResult analysePair(DependencePair*);

        void getMemRefInstrs(const Loop *, SmallVectorImpl<Instruction*> &);
        Value *getPointerOperand(Value *);
        const SCEV *getZeroSCEV(ScalarEvolution *);
        bool isMemRefInstr(const Value *);
        bool isLoadOrStoreInst(Value *);
        bool isDependencePair(const Value *A, const Value *B);
        bool depends(Value *, Value *);
        AliasAnalysis::AliasResult underlyingObjectsAlias(AliasAnalysis *, const Value *, const Value *);
        void PrintLoopInfo(raw_ostream &OS, CS551A2 *, const Loop*);

        void partition();

        // END COPY-PASTE }}}

        static char ID;
    public:
        CS551A2() : LoopPass(ID) {}
        virtual bool runOnLoop(Loop *LP, LPPassManager &LPM);
        virtual void print(raw_ostream &O, const Module *M) const;
        virtual void releaseMemory();
        virtual void getAnalysisUsage(AnalysisUsage&) const;

    private:
        FoldingSet<DependencePair> Pairs;
        BumpPtrAllocator PairAllocator;
        AliasAnalysis *AA;
        /*!
         @brief Contains the Scalar Evaluator
         */
        ScalarEvolution *SE;
        LoopInfo *LI;
        Loop *L;
    };
}

// I do not understand what this does;
// it is just copy-and-pasted from the example :-(
char CS551A2::ID = 0;
static RegisterPass<CS551A2> X("cs551a2", "Loop Dependency Magick");

void CS551A2::partition()
{
    // S := set(subscript pairs)
    // P[out] := set(separable or minimal coupled)
    // Np := len(P)
}

void CS551A2::print(raw_ostream &O, const Module *M) const
{
    O << "hello, I am a CS551a2 module\n";
    // we have to throw away the const qualifier in order
    // to get back into the function graph where everything
    // is NOT declared "const". We have no control over the
    // signature of _this_ function because it is coming from LoopPass
    CS551A2 *nonConst = (CS551A2*)this;
    nonConst->PrintLoopInfo(O, nonConst, this->L);
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/*
INITIALIZE_PASS_BEGIN(LoopDependenceAnalysis, "lda",
                "Loop Dependence Analysis", false, true)
INITIALIZE_PASS_DEPENDENCY(ScalarEvolution)
INITIALIZE_AG_DEPENDENCY(AliasAnalysis)
INITIALIZE_PASS_END(LoopDependenceAnalysis, "lda",
                "Loop Dependence Analysis", false, true)
char LoopDependenceAnalysis::ID = 0;
*/

//===----------------------------------------------------------------------===//
//                             Utility Functions
//===----------------------------------------------------------------------===//

bool CS551A2::isMemRefInstr(const Value *V) {
  const Instruction *I = dyn_cast<const Instruction>(V);
  return I && (I->mayReadFromMemory() || I->mayWriteToMemory());
}

void CS551A2::getMemRefInstrs(const Loop *L, SmallVectorImpl<Instruction*> &Memrefs) {
  for (Loop::block_iterator b = L->block_begin(), be = L->block_end();
       b != be; ++b)
    for (BasicBlock::iterator i = (*b)->begin(), ie = (*b)->end();
         i != ie; ++i)
      if (isMemRefInstr(i))
        Memrefs.push_back(i);
}

/*! @brief Decides if the presented Value is a Load or Store instruction.
 */
bool CS551A2::isLoadOrStoreInst(Value *I) {
  // Returns true if the load or store can be analyzed. Atomic and volatile
  // operations have properties which this analysis does not understand.
  if (LoadInst *LI = dyn_cast<LoadInst>(I))
    return LI->isUnordered();
  else if (StoreInst *SI = dyn_cast<StoreInst>(I))
    return SI->isUnordered();
  return false;
}

/*! @brief Returns the pointer operand of the given Load or Store instruction.
 If you provide an instruction that is neither a Load or Store, you'll be sorry.
 */
Value *
CS551A2::getPointerOperand(Value *I) {
  if (LoadInst *i = dyn_cast<LoadInst>(I))
    return i->getPointerOperand();
  if (StoreInst *i = dyn_cast<StoreInst>(I))
    return i->getPointerOperand();
  llvm_unreachable("Value is no load or store instruction!");
  // Never reached.
  return 0;
}

AliasAnalysis::AliasResult
CS551A2::underlyingObjectsAlias(AliasAnalysis *AA, const Value *A, const Value *B) {
  const Value *aObj = GetUnderlyingObject(A);
  const Value *bObj = GetUnderlyingObject(B);
  return AA->alias(aObj, AA->getTypeStoreSize(aObj->getType()),
                   bObj, AA->getTypeStoreSize(bObj->getType()));
}

const SCEV *
CS551A2::getZeroSCEV(ScalarEvolution *SE) {
  return SE->getConstant(Type::getInt32Ty(SE->getContext()), 0L);
}

//===----------------------------------------------------------------------===//
//                             Dependence Testing
//===----------------------------------------------------------------------===//

bool CS551A2::findOrInsertDependencePair(Value *A, Value *B, DependencePair *&P) {
  void *insertPos = 0;
  FoldingSetNodeID id;
  id.AddPointer(A);
  id.AddPointer(B);

  P = this->Pairs.FindNodeOrInsertPos(id, insertPos);
  if (P) return true;

  P = new (PairAllocator)DependencePair(id, A, B);
  this->Pairs.InsertNode(P, insertPos);
  return false;
}

/*! @brief Gets the Loops where "S" is a Loop invariant,
 returning the answer in "Loops".
 */
void CS551A2::getLoops(const SCEV *S, DenseSet<const Loop*>* Loops) const {
  // Refactor this into an SCEVVisitor, if efficiency becomes a concern.
    for (const Loop *L = this->L; L != 0; L = L->getParentLoop()) {
        if (!SE->isLoopInvariant(S, L)) {
            Loops->insert(L);
        }
    }
}

bool CS551A2::isLoopInvariant(const SCEV *S) const {
  DenseSet<const Loop*> loops;
  getLoops(S, &loops);
  return loops.empty();
}

bool CS551A2::isAffine(const SCEV *S) const {
    const SCEVAddRecExpr *rec = dyn_cast<SCEVAddRecExpr>(S);
    bool invar = isLoopInvariant(S);
    DEBUG(dbgs() << "SCEV[" << *S << "]:isLoopInvariant? " << invar << "\n");
    bool isAddRecExpr = NULL != rec;
    DEBUG(dbgs() << "SCEV[" << *S << "]:isAndRecExpr? " << isAddRecExpr << "\n");
    bool isAffine = (rec && rec->isAffine());
    DEBUG(dbgs() << "SCEV[" << *S << "]:is SE Affine? " << isAffine << "\n");
    bool result = (invar || isAffine);
    DEBUG(dbgs() << "::isAffine(" << *S << ") <- " << result << "\n");
    return result;
}

CS551A2::Affine*
CS551A2::getAffine(const Value *V)
{
    /// use a tmp because we are going to reassign it a lot
    const Value *tmpV = V;
    Affine* result = new Affine;
    result->coefficient = 1;
    result->constant = 0;
    while (NULL != tmpV) {
        const Instruction *bI = dyn_cast<const Instruction>(tmpV);
        if (! bI) {
            dbgs() << "getAffine bottomed out on " << *tmpV << "\n";
            break;
        }
        dbgs() << "switching on " << *bI << "\n";

        if (const SExtInst *sext = dyn_cast<const SExtInst>(tmpV)) {
            dbgs() << "sign-ext!\n";
            tmpV = sext->getOperand(0);
        } else if (const AddOperator *addOp = dyn_cast<const AddOperator>(tmpV)) {
            const Value *op1 = addOp->getOperand(0);
            const Value *op2 = addOp->getOperand(1);
            dbgs() << "add! OP("<< *op1 << "),OP(" << *op2 << ")\n";

            /// kill the iteration unless we figure something out
            tmpV = NULL;
            bool nextOp1 = false;
            bool nextOp2 = false;
            if (const ConstantInt *ci = dyn_cast<const ConstantInt>(op1)) {
                dbgs() << "op1 is constant: " << *ci << "\n";
                result->constant = ci->getLimitedValue();
            } else if (const MulOperator *mulOp = dyn_cast<const MulOperator>(op1)) {
                dbgs() << "op1 is a mult" << *mulOp <<"\n";
                nextOp1 = true;
            } else if (const LoadInst *li = dyn_cast<const LoadInst>(op1)) {
                dbgs() << "op1 is loaded" << *li <<"\n";
                nextOp1 = true;
            }
            if (const ConstantInt *ci = dyn_cast<const ConstantInt>(op2)) {
                dbgs() << "op2 is constant: " << *ci << "\n";
                result->constant = ci->getLimitedValue();
            } else if (const MulOperator *mulOp = dyn_cast<const MulOperator>(op2)) {
                dbgs() << "op2 is a mult" << *mulOp <<"\n";
                nextOp2 = true;
            } else if (const LoadInst *li = dyn_cast<const LoadInst>(op2)) {
                dbgs() << "op2 is loaded" << *li <<"\n";
                nextOp2 = true;
            }
            if (nextOp1 && !nextOp2) {
                dbgs() << "Looks like we'll proceed with op1\n";
                tmpV = op1;
            } else if (!nextOp1 && nextOp2) {
                dbgs() << "Looks like we'll proceed with op2\n";
                tmpV = op2;
            }
        } else if (const MulOperator *mulOp = dyn_cast<const MulOperator>(tmpV)) {
            const Value *op1 = mulOp->getOperand(0);
            const Value *op2 = mulOp->getOperand(1);
            dbgs() << "multiply! OP("<< *op1 << "),OP(" << *op2 << ")\n";

            /// my kingdom for a lambda here :-(

            /// kill the iteration unless we figure something out
            tmpV = NULL;
            bool nextOp1 = false;
            bool nextOp2 = false;
            if (const ConstantInt *ci = dyn_cast<const ConstantInt>(op1)) {
                dbgs() << "op1 is constant: " << *ci << "\n";
                result->coefficient = ci->getLimitedValue();
            } else if (const AddOperator *addOp = dyn_cast<const AddOperator>(op1)) {
                dbgs() << "op1 is an add" << *addOp <<"\n";
                nextOp1 = true;
            } else if (const LoadInst *li = dyn_cast<const LoadInst>(op1)) {
                dbgs() << "op1 is loaded" << *li <<"\n";
                nextOp1 = true;
            }
            if (const ConstantInt *ci = dyn_cast<const ConstantInt>(op2)) {
                dbgs() << "op2 is constant: " << *ci << "\n";
                result->coefficient = ci->getLimitedValue();
            } else if (const AddOperator *addOp = dyn_cast<const AddOperator>(op2)) {
                dbgs() << "op2 is an add" << *addOp <<"\n";
                nextOp2 = true;
            } else if (const LoadInst *li = dyn_cast<const LoadInst>(op2)) {
                dbgs() << "op2 is loaded" << *li <<"\n";
                nextOp2 = true;
            }
            if (nextOp1 && !nextOp2) {
                dbgs() << "Looks like we'll proceed with op1\n";
                tmpV = op1;
            } else if (!nextOp1 && nextOp2) {
                dbgs() << "Looks like we'll proceed with op2\n";
                tmpV = op2;
            }
        } else if (const LoadInst *li = dyn_cast<const LoadInst>(tmpV)) {
            const Value *theVar = li->getOperand(0);
            dbgs() << "Victory, we hit a memory inst with " << *li
                << " with variable: "<< *theVar << "\n\n";
            result->index = theVar;
            tmpV = NULL;
        } else {
            dbgs() << "Sorry, I don't know what to make of " << *tmpV << "\n";
            tmpV = NULL;
        }
    }
    if (result->index) {
        return result;
    } else {
        /// we didn't find an index variable, so not affine
        return NULL;
    }
}

bool CS551A2::isZIVPair(const SCEV *A, const SCEV *B) const {
    DEBUG(dbgs() << "isZIVPair(" << *A << "," << *B << ") ...\n");
    bool result;
    result = isLoopInvariant(A) && isLoopInvariant(B);
    DEBUG(dbgs() << "isZIVPair(" << *A << "," << *B << ") <- " << result << "\n");
    return result;
}

bool CS551A2::isSIVPair(const SCEV *A, const SCEV *B) const {
  DenseSet<const Loop*> loops;
  getLoops(A, &loops);
  getLoops(B, &loops);
  return loops.size() == 1;
}

CS551A2::DependenceResult
CS551A2::analyseZIV(const SCEV *A, const SCEV *B) const {
  assert(isZIVPair(A, B) && "Attempted to ZIV-test non-ZIV SCEVs!");
  return A == B ? Dependent : Independent;
}

CS551A2::DependenceResult
CS551A2::analyseSIV(const SCEV *A, const SCEV *B) const {
  return Unknown; // TODO: Implement.
}

CS551A2::DependenceResult
CS551A2::analyseMIV(const SCEV *A, const SCEV *B) const {
  return Unknown; // TODO: Implement.
}

CS551A2::DependenceResult
CS551A2::analyseSubscript(const SCEV *A, const SCEV *B) const {
    DEBUG(dbgs() << "\n\nTesting subscript: A(" << *A << "), B(" << *B << ")\n");

  if (A == B) {
    DEBUG(dbgs() << "  -> [D] same SCEV\n");
    return Dependent;
  }

    if (!isAffine(A)) {
        DEBUG(dbgs() << "  -> [?] A is not affine\n");
        return Unknown;
    }
    if (!isAffine(B)) {
        DEBUG(dbgs() << "  -> [?] B is not affine\n");
        return Unknown;
    }

  if (isZIVPair(A, B))
    return analyseZIV(A, B);

  if (isSIVPair(A, B))
    return analyseSIV(A, B);

  return analyseMIV(A, B);
}

static bool isStrongSIV(CS551A2::Subscript *S)
{
    if (S->A->index != S->B->index) {
        dbgs() << "Not Strong due to differing indices\n";
        return false;
    }
    bool result = S->A->coefficient == S->B->coefficient;
    return result;
}

CS551A2::DependenceResult
CS551A2::analysePair(DependencePair *P) {
  DEBUG(dbgs() << "Analysing:\n" << *P->A << "\n" << *P->B << "\n");

  // We only analyse loads and stores but no possible memory accesses by e.g.
  // free, call, or invoke instructions.
  if (!isLoadOrStoreInst(P->A) || !isLoadOrStoreInst(P->B)) {
    DEBUG(dbgs() << "--> [?] no load/store\n");
    return Unknown;
  }

  Value *aPtr = getPointerOperand(P->A);
  Value *bPtr = getPointerOperand(P->B);

  switch (underlyingObjectsAlias(AA, aPtr, bPtr)) {
  case AliasAnalysis::MayAlias:
  case AliasAnalysis::PartialAlias:
    // We can not analyse objects if we do not know about their aliasing.
    DEBUG(dbgs() << "---> [?] may alias\n");
    return Unknown;

  case AliasAnalysis::NoAlias:
    // If the objects noalias, they are distinct, accesses are independent.
    DEBUG(dbgs() << "---> [I] no alias\n");
    return Independent;

  case AliasAnalysis::MustAlias:
    break; // The underlying objects alias, test accesses for dependence.
  }

  /// GEP is a Get Element Ptr instruction
  const GEPOperator *aGEP = dyn_cast<GEPOperator>(aPtr);
  const GEPOperator *bGEP = dyn_cast<GEPOperator>(bPtr);

  if (!aGEP || !bGEP)
    return Unknown;

    DEBUG(dbgs() << "GEP-A := " << *aGEP << "\n" << "GEP-B := " << *bGEP << "\n");

  // FIXME: Is filtering coupled subscripts necessary?

    // Collect GEP operand pairs (FIXME: use GetGEPOperands from BasicAA), adding
    // trailing zeroes to the smaller GEP, if needed.

    const SCEV *SCEV_ZERO = getZeroSCEV(SE);

    for(GEPOperator::const_op_iterator aIdx = aGEP->idx_begin(),
        aEnd = aGEP->idx_end(),
        bIdx = bGEP->idx_begin(),
        bEnd = bGEP->idx_end();
        aIdx != aEnd && bIdx != bEnd;
        aIdx += (aIdx != aEnd), bIdx += (bIdx != bEnd))
    {
        const Value *aValue = *aIdx;
        const Value *bValue = *bIdx;

        DEBUG(dbgs() << "aIdx := " << *aValue << "\tType:" << *(aValue->getType()) << "\n");
        const SCEV* aSCEV = (aIdx != aEnd)
            ? SE->getSCEV(*aIdx)
            : SCEV_ZERO;
        
        DEBUG(dbgs() << "bIdx := " << *bValue << "\tType:" << *(bValue->getType()) << "\n");
        const SCEV* bSCEV = (bIdx != bEnd)
            ? SE->getSCEV(*bIdx)
            : SCEV_ZERO ;
        if (aSCEV == SCEV_ZERO && bSCEV == SCEV_ZERO) {
            DEBUG(dbgs() << "skipping this round because both SCEV are ZERO\n");
            continue;
        }
        DEBUG(dbgs() << "SCEV-A := " << *aSCEV << "\n" << "SCEV-B := " << *bSCEV << "\n");
        if (SE->isLoopInvariant(aSCEV, this->L)) {
            dbgs() << "Hey, aSCEV is loop-invariant\n";
        }
        if (SE->isLoopInvariant(bSCEV, this->L)) {
            dbgs() << "Hey, bSCEV is loop-invariant\n";
        }

        const Affine *aAffine;
        const Affine *bAffine;
        if (aSCEV != SCEV_ZERO) {
            if ((aAffine = getAffine(aValue))) {
                dbgs() << "Successfully affined aValue := "
                << aAffine->coefficient
                << aAffine->index->getName()
                << "+" << aAffine->constant << "\n";
            }
        }
        if (bSCEV != SCEV_ZERO) {
            if ((bAffine = getAffine(bValue))) {
                dbgs() << "Successfully affined bValue := "
                << bAffine->coefficient
                << bAffine->index->getName()
                << "+" << bAffine->constant << "\n";
            }
        }

        Subscript *sub;
        if (aAffine && bAffine) {
            sub = new Subscript;
            sub->A = aAffine;
            sub->B = bAffine;
            if (sub->A->index == sub->B->index) {
                dbgs() << "Delta(" << sub->A->index->getName()
                    << ")="
                    << (sub->B->constant - sub->A->constant)
                    << "\n";
                dbgs() << "Strong? " << ::isStrongSIV(sub) << "\n";
                return Dependent;
            } else {
                return Independent;
            }
        }
    }

    // We successfully analysed all subscripts but failed to prove independence.
    return Dependent;
}

bool CS551A2::depends(Value *A, Value *B) {
  assert(isDependencePair(A, B) && "Values form no dependence pair!");

  DependencePair *p;
  if (!findOrInsertDependencePair(A, B, p)) {
    switch (p->Result = analysePair(p)) {
    case Dependent:
            errs() << "DPair("<< *A <<","<< *B <<") is Dependent\n";
            break;
    case Independent:
            errs() << "DPair("<< *A <<","<< *B <<") is IN-dependent\n";
            break;
    case Unknown:
            errs() << "DPair("<< *A <<","<< *B <<") is UNKNOWN\n";
            break;
    }
  }
  return p->Result != Independent;
}

//===----------------------------------------------------------------------===//
//                   LoopDependenceAnalysis Implementation
//===----------------------------------------------------------------------===//

bool CS551A2::runOnLoop(Loop *L, LPPassManager &) {
  this->L = L;
  AA = &getAnalysis<AliasAnalysis>();
  SE = &getAnalysis<ScalarEvolution>();
  DEBUG(this->dump());
  return false;
}

void CS551A2::releaseMemory() {
  Pairs.clear();
  PairAllocator.Reset();
}

void CS551A2::getAnalysisUsage(AnalysisUsage &AU) const {
  AU.setPreservesAll();
  AU.addRequiredTransitive<AliasAnalysis>();
  AU.addRequiredTransitive<ScalarEvolution>();
}

bool CS551A2::isDependencePair(const Value *A, const Value *B) {
    return
        isMemRefInstr(A) &&
        isMemRefInstr(B) &&
    (cast<const Instruction>(A)->mayWriteToMemory() ||
     cast<const Instruction>(B)->mayWriteToMemory());
}

void CS551A2::PrintLoopInfo(raw_ostream &OS, CS551A2 *self, const Loop* LP) {
  if (!LP || !LP->empty()) return; // ignore non-innermost loops

  SmallVector<Instruction*, 8> memrefs;
  getMemRefInstrs(LP, memrefs);

  OS << "Loop at depth " << LP->getLoopDepth() << ", header block: ";
  WriteAsOperand(OS, LP->getHeader(), false);
  OS << "\n";

  OS << "  Load/store instructions: " << memrefs.size() << "\n";
  for (SmallVector<Instruction*, 8>::const_iterator x = memrefs.begin(),
       end = memrefs.end(); x != end; ++x)
    OS << "\t" << (x - memrefs.begin()) << ": " << **x << "\n";

  OS << "  Pairwise dependence results:\n";
  for (SmallVector<Instruction*, 8>::const_iterator x = memrefs.begin(),
       end = memrefs.end(); x != end; ++x)
    for (SmallVector<Instruction*, 8>::const_iterator y = x + 1;
         y != end; ++y)
      if (isDependencePair(*x, *y))
        OS << "\t" << (x - memrefs.begin()) << "," << (y - memrefs.begin())
           << ": " << (depends(*x, *y) ? "dependent" : "independent")
           << "\n";
}
