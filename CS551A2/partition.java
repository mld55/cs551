import java.util.ArrayList;
import java.util.List;

public class partition
{
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
        public Subscript(IndexEntry lhs, IndexEntry rhs) {
            this.lhs = lhs;
            this.rhs = rhs;
        }
        public boolean hasIndex(String idx) {
            return lhs.hasIndex(idx) || rhs.hasIndex(idx);
        }
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
        // A[I+1][I][K] = A[I][J][K]
        final String[] indexes = { "I","J","K" };
        final List<List<Subscript>> partitions
            = new ArrayList<List<Subscript>>();
        Subscript s0 = new Subscript(
new IndexEntry(new int[] { 1, 0, 0 }, indexes, 1),
new IndexEntry(new int[] { 1, 0, 0 }, indexes, 0));
        Subscript s1 = new Subscript(
new IndexEntry(new int[] { 1, 0, 0 }, indexes, 0),
new IndexEntry(new int[] { 0, 1, 0 }, indexes, 0));
        Subscript s2 = new Subscript(
new IndexEntry(new int[] { 0, 0, 1 }, indexes, 0),
new IndexEntry(new int[] { 0, 0, 1 }, indexes, 0));
        final List<Subscript> subs
            = new ArrayList<Subscript>(3);
        subs.add( s0 );
        subs.add( s1 );
        subs.add( s2 );
        System.out.println("Subs := "+subs);
        partition( subs, partitions, indexes );
        System.out.println("PART := "+partitions);
    }
}
