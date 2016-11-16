package au.edu.wehi.idsv.alignment;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.UUID;

import com.google.common.io.Files;

import au.edu.wehi.idsv.Defaults;
import htsjdk.samtools.util.Log;

public class AlignerFactory {
	private static final Log log = Log.getInstance(AlignerFactory.class);
	private static final String SSW_JNI_JAR_LOCATION = "/libsswjni.so";
	private static boolean sswjniLoaded;
	private static final Aligner defaultAligner;
    static {
    	if (!Defaults.NO_LIBSSW) {
    		try {
            	System.loadLibrary("sswjni");
            	sswjniLoaded = true;
            } catch (UnsatisfiedLinkError e) {
            	sswjniLoaded = false;
            	log.debug("Unable to load sswjni library");
            }
            try {
            	if (!sswjniLoaded) {
            		File tmp = new File(System.getProperty("java.io.tmpdir"), UUID.randomUUID().toString() + "-libsswjni.so");
                	File destination = new File(System.getProperty("java.io.tmpdir"), "libsswjni.so");
                	try {
                		unpacksswjni(tmp);
                		Files.move(tmp, destination);
                		System.load(destination.getAbsolutePath());
                		sswjniLoaded = true;
                	} catch (IOException e) {
                		log.error("Unable to extract sswjni native library to " + destination.toString());
                	}
                }
            } catch (UnsatisfiedLinkError e) {
            	sswjniLoaded = false;
            }
    	} else {
    		sswjniLoaded = false;
    	}
        defaultAligner = create(1, -4, -4, 6, 1); // bwa mem
        // defaultAligner = create(2, -6, -1, 5, 3); // bowtie2
        if (sswjniLoaded) {
        	try {
        		log.debug("Testing JNI alignment");
	        	defaultAligner.align_smith_waterman(new byte[] { 'A' } , new byte[] { 'A' });
	        	log.info("sswjni library loading successful.");
        	} catch (Exception e) {
        		log.error(e, "sswjni JNI call failure. Attempting to recover.");
        		sswjniLoaded = false;
        	}
        }
        if (!sswjniLoaded) {
        	log.error("Unable to use sswjni library - assembly will be very slow. Please ensure libsswjni for your OS and architecture can be found on java.library.path");
        }
    }
    private static void unpacksswjni(File destination) throws IOException {
    	if (destination.exists() && destination.length() == AlignerFactory.class.getResource(SSW_JNI_JAR_LOCATION).getFile().length()) {
    		log.debug("Found " + destination.toString());
    		return;
    	}
        try (InputStream reader = AlignerFactory.class.getResourceAsStream(SSW_JNI_JAR_LOCATION)) {
            try (FileOutputStream writer = new FileOutputStream(destination)) {
                byte[] buffer = new byte[8192];
                int bytesRead = 0;
                while ((bytesRead = reader.read(buffer)) != -1) {
                    writer.write(buffer, 0, bytesRead);
                }
            }
        }
        log.debug("Created " + destination.toString());
    }
	public static Aligner create(int match, int mismatch, int ambiguous, int gapOpen, int gapExtend) {
		if (sswjniLoaded) {
			return new SswJniAligner(match, mismatch, ambiguous, gapOpen, gapExtend);
		} else {
			return new JAlignerAligner(match, mismatch, ambiguous, gapOpen, gapExtend);
		}
	}
	public static Aligner create() {
		return defaultAligner;
	}
}
