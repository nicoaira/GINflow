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
}
