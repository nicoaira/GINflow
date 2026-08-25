/**
 * Parse Nextflow CLI/config booleans.
 *
 * Nextflow 26 (strict parser) keeps {@code --flag false} as the string
 * {@code "false"}, which is Groovy-truthy. Reassigning {@code params.flag}
 * after that is ignored. Call this at every boolean decision site.
 */
class BooleanParam {

    static boolean parse(value) {
        return parse(value, false)
    }

    static boolean parse(value, boolean defaultValue) {
        if (value == null) {
            return defaultValue
        }
        if (value instanceof Boolean) {
            return (boolean) value
        }
        def text = value.toString().trim().toLowerCase()
        if (text in ['true', '1', 'yes', 'y', 'on']) {
            return true
        }
        if (text in ['false', '0', 'no', 'n', 'off', '']) {
            return false
        }
        throw new IllegalArgumentException("must be true or false, got '${value}'")
    }

    static boolean rerankEnabled(exactRerank, hnswlibRerank) {
        return parse(exactRerank, true) || parse(hnswlibRerank, false)
    }

    static boolean codeQuantized(quantize) {
        def text = quantize == null ? 'none' : quantize.toString().trim().toLowerCase()
        return text in ['pq', 'opq']
    }

    /**
     * Cosine {@code --seed_min_similarity} applies to original-window scores.
     * PQ/OPQ ADC is not cosine (on the 30k Rfam index the best ADC was ~0.55),
     * so un-reranked code search must not use that cutoff.
     */
    static String annMinSimilarity(exactRerank, hnswlibRerank, quantize, seedMinSimilarity) {
        if (rerankEnabled(exactRerank, hnswlibRerank) || codeQuantized(quantize)) {
            return '-inf'
        }
        return String.valueOf(seedMinSimilarity)
    }
}
