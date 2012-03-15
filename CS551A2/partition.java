import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class partition
{
    enum DependencyDirection
    {
        LT("<"), EQ("="), GT(">"), STAR("*");
        private DependencyDirection(String str) {
            this.str = str;
        }
        @Override
        public String toString() {
            return str;
        }
        private String str;
    }
    static class IndexEntry
    {
        public IndexEntry(int[] co, String[] indexes, int constant) {
            if (co != null && indexes != null
                    && co.length != indexes.length) {
                throw new IllegalArgumentException(
"You must provide coefficients for all your indexes, even if 0 or 1");
            }
            this.coefficients = co;
            this.indexes = indexes;
            this.constant = constant;
        }
        /**
         * Contains the coefficient before the index variable, or 1 if
         * does not have one.  Please note this can be zero, if an
         * index variable does not appear in this index entry.
        */
        int[] coefficients;
        /** Contains the name of the index variable. */
        String[] indexes;
        /** Contains any constant for this index entry. */
        int constant;
        /**
         * Returns true if this index entry has a non-zero coefficient index variable you provide.
         */
        public boolean hasIndex(String idx) {
            if (null == indexes || 0 == indexes.length) {
                return false;
            }
            for (int i = 0; i < indexes.length; i++) {
                final String index = indexes[i];
                final int coeff = getCoefficient(i);
                if (0 == coeff) {
                    continue;
                }
                if (idx.equals(index)) {
                    return true;
                }
            }
            return false;
        }
        /** Returns the coefficient for index at <tt>pos</tt>, or
            1 if the user did not specify coefficients. */
        private int getCoefficient(int pos) {
            final int result;
            if (null == coefficients) {
                result = 1;
            } else {
                if (pos >= coefficients.length && pos >= indexes.length) {
                    throw new IllegalArgumentException("There is no such coefficient");
                }
                result = coefficients[ pos ];
            }
            return result;
        }
        @Override public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("Index(");
            if (null != indexes) {
                // controlling where to place the formatting comma because
                // the "i" variable is not enough information given that some
                // indexes can be zero-ed out
                boolean first = true;
                for (int i = 0; i < indexes.length; i++) {
                    final int coeff = getCoefficient(i);
                    if (0 == coeff) {
                        // this index is zeroed
                        continue;
                    }
                    // wait until we have tested the coeff before formatting
                    if (! first) {
                        sb.append(",");
                    }
                    if (1 == coeff) {
                        // don't print anything, otherwise it'll be "1I+1"
                        // which looks horrible
                    } else {
                        sb.append(coeff);
                    }
                    sb.append( indexes[i] );
                    first = false;
                }
            }
            if (0 != constant) {
                sb.append( constant < 0 ? "-" : "+" )
                    .append( constant );
            }
            sb.append(")");
            return sb.toString();
        }
    }

    static class Subscript
    {
        public Subscript(String[] indexes, IndexEntry lhs, IndexEntry rhs) {
            if (null != lhs.indexes && null != rhs.indexes) {
                if (lhs.indexes.length != rhs.indexes.length) {
                    throw new IllegalArgumentException(
                            "Index Entry indexes must have equal lengths");
                }
            }
            this.indexes = indexes;
            this.lhs = lhs;
            this.rhs = rhs;
        }

        public boolean hasIndex(String idx) {
            return lhs.hasIndex(idx) || rhs.hasIndex(idx);
        }

        public int getIndexCount() {
            int count = 0;
            for (String i : indexes) {
                if (lhs.hasIndex(i) || rhs.hasIndex(i)) {
                    count++;
                }
            }
            return count;
        }

        public boolean isZIV() {
            return 0 == getIndexCount();
        }

        public boolean isSIV() {
            return 1 == getIndexCount();
        }

        public boolean isMIV() {
            return 1 < getIndexCount();
        }

        public boolean isStrongSIV() {
            if (! isSIV() ) return false;
            for (int i = 0; i < indexes.length; i++) {
                final String idx = indexes[i];
                if (lhs.hasIndex(idx)) {
                    return lhs.getCoefficient(i) ==
                        rhs.getCoefficient(i);
                }
            }
            return false;
        }

        public boolean isWeakSIV() {
            if (! isSIV() ) return false;
            for (int i = 0; i < indexes.length; i++) {
                final String idx = indexes[i];
                // the fact that it *has* the co-eff implies it is non-zero
                if (lhs.hasIndex(idx)) {
                    return lhs.getCoefficient(i) != rhs.getCoefficient(i);
                }
            }
            return false;
        }

        public boolean isWeakZeroSIV() {
            if (! isSIV() ) return false;
            for (int i = 0; i < indexes.length; i++) {
                final String idx = indexes[i];
                if (lhs.hasIndex(idx)) {
                    return !rhs.hasIndex(idx);
                } else if (rhs.hasIndex(idx)) {
                    return !lhs.hasIndex(idx);
                }
            }
            return false;
        }

        public List<DependencyDirection> getDirectionVector() {
            final List<DependencyDirection> results
                = new ArrayList<DependencyDirection>(indexes.length);
            for (int i = 0; i < indexes.length; i++) {
                final String idx = indexes[i];
                if (!( lhs.hasIndex(idx) && rhs.hasIndex(idx))) {
                    results.add( DependencyDirection.STAR );
                    continue;
                }
                int l_coeff = lhs.getCoefficient(i);
                int r_coeff = rhs.getCoefficient(i);
                DependencyDirection dd;
                if (l_coeff < r_coeff) {
                    dd = DependencyDirection.GT;
                } else if (r_coeff > r_coeff) {
                    dd = DependencyDirection.LT;
                } else {
                    // co-eff is equal
                    if (lhs.constant < rhs.constant) {
                        dd = DependencyDirection.GT;
                    } else if (lhs.constant > rhs.constant) {
                        dd = DependencyDirection.LT;
                    } else { // if (lhs.constant == rhs.constant) {
                        dd = DependencyDirection.EQ;
                    }
                }
                results.add( dd );
            }
            return results;
        }

        String[] indexes;
        IndexEntry lhs;
        IndexEntry rhs;
        @Override public String toString() {
            return String.format("Subscript<%s, %s>", lhs, rhs);
        }
    }

    /**
     * Partitions the provided Subscripts into coupled and non-coupled lists.
     * @param subs the list of subscripts to partition
     * @param partitions the <u>output</u> list of partitions, meaning
     * you should pass in a fresh List and I will fill it for you.
     * @param indexVars the list of outer-most to inner-most index
     * variable names, which should match up with the index variables
     * found in the <tt>Subscript</tt> instances
     */
    public static void partition( List<Subscript> subs, List<List<Subscript>> partitions, String[] indexVars)
    {
        for (final Subscript s : subs) {
            final List<Subscript> list = new ArrayList<Subscript>(1);
            list.add( s );
            partitions.add( list );
        }
        for (final String idx : indexVars) {
            int k = -1;
            // we are using "downto" because if we kill
            // an entry, we can continue the for loop,
            // which is not true when going "up"
            for (int j = partitions.size() - 1; j >= 0; j--) {
                final List<Subscript> list = partitions.get(j);
                for (int x = list.size() - 1; x >= 0; x--) {
                    final Subscript it = list.get(x);
                    if (it.hasIndex(idx)) {
                        if (-1 == k) {
                            k = j;
                        } else {
                            list.remove(it);
                            partitions.get(k).add(it);
                        }
                    }
                }
                if (list.isEmpty()) {
                    partitions.remove(list);
                }
            }
        }
    }

    public static void main(String[] args) throws Exception {
        final List<List<Subscript>> partitions
            = new ArrayList<List<Subscript>>();
/*
        // A[I+1][I][K][5] = A[I][J][K][8]
        final String[] indexes = { "I","J","K" };
        Subscript s0 = new Subscript( indexes,
            new IndexEntry(new int[] { 2, 0, 0 }, indexes, 1),
            new IndexEntry(new int[] { 0, 0, 0 }, indexes, 1));
        Subscript s1 = new Subscript( indexes,
            new IndexEntry(new int[] { 1, 0, 0 }, indexes, 0),
            new IndexEntry(new int[] { 0, 1, 0 }, indexes, 0));
        Subscript s2 = new Subscript( indexes,
            new IndexEntry(new int[] { 0, 0, 1 }, indexes, 0),
            new IndexEntry(new int[] { 0, 0, 1 }, indexes, 0));
        Subscript s3 = new Subscript( indexes,
            new IndexEntry(new int[] { 0, 0, 0 }, indexes, 5),
            new IndexEntry(new int[] { 0, 0, 0 }, indexes, 8));
*/
/*
        // A[i+1][i] = A[i][i+1]
        final String[] indexes = { "I" };
        Subscript s0 = new Subscript( indexes,
                new IndexEntry(new int[] { 1 }, indexes, 1),
                new IndexEntry(new int[] { 1 }, indexes, 0));
        Subscript s1 = new Subscript( indexes,
                new IndexEntry(new int[] { 1 }, indexes, 0),
                new IndexEntry(new int[] { 1 }, indexes, 1));
*/
        // A[i+1][i+2] = A[i][i]
        final String[] indexes = { "I" };
        Subscript s0 = new Subscript( indexes,
                new IndexEntry(new int[] { 1 }, indexes, 1),
                new IndexEntry(new int[] { 1 }, indexes, 0));
        Subscript s1 = new Subscript( indexes,
                new IndexEntry(new int[] { 1 }, indexes, 2),
                new IndexEntry(new int[] { 1 }, indexes, 0));

        final List<Subscript> subs = Arrays.asList( s0, s1 );
        System.out.println("Subs := "+subs);
        partition( subs, partitions, indexes );
        System.out.println("PART := "+partitions);
        for (final List<Subscript> part : partitions) {
            System.out.println("--PART");
            for (final Subscript sub : part) {
                System.out.println("SUB="+sub);
                System.out.println("\tZIV? "+sub.isZIV());
                System.out.println("\tSIV? "+sub.isSIV());
                System.out.println("\t\tStrong-SIV? "+sub.isStrongSIV());
                System.out.println("\t\tWeak-SIV? "+sub.isWeakSIV());
                System.out.println("\t\tWeakZero-SIV? "+sub.isWeakZeroSIV());
                System.out.println("\tMIV? "+sub.isMIV());
                System.out.println("\tDIR="+sub.getDirectionVector());
            }
        }
    }
}