import groovy.json.JsonSlurper

class PannotatorUtils {
    public static String get_bakta_db_type(String bakta_db) {
        def version_json_path = new File(bakta_db, "version.json")

        if (!version_json_path.exists()) {
            return null
        }

        try {
            def version_json = new JsonSlurper().parse(version_json_path)
            return version_json?.type as String
        } catch (Exception e) {
            return null
        }
    }
}
