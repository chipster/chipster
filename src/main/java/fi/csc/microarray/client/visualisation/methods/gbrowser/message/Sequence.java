package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

/**
 * Utility class for DNA sequence operations.
 * 
 */
public class Sequence {
    
    /**
     * Find a complement for given nucleotide sequence.
     * 
     * @param seq
     * @return
     */
    public static String complement(String seq) {

        StringBuffer buf = new StringBuffer(seq);

        for (int j = 0; j < seq.length(); j++) {
            switch (buf.charAt(j)) {
            case 'A':
                buf.setCharAt(j, 'T');
                break;
            case 'C':
                buf.setCharAt(j, 'G');
                break;
            case 'G':
                buf.setCharAt(j, 'C');
                break;
            case 'T':
                buf.setCharAt(j, 'A');
                break;
            case 'a':
                buf.setCharAt(j, 't');
                break;
            case 'c':
                buf.setCharAt(j, 'g');
                break;
            case 'g':
                buf.setCharAt(j, 'c');
                break;
            case 't':
                buf.setCharAt(j, 'a');
                break;
            }
        }

        return buf.toString();
    }
    
}
