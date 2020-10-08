package au.edu.wehi.idsv.util;

public abstract class FilenameUtil {
    public static String stripInvalidFilenameCharacters(String filename) {
        return filename.replaceAll("[|/\\<>?*\"\\[\\],:;]", "");
    }
}
