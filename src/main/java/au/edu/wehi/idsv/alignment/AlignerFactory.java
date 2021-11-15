package au.edu.wehi.idsv.alignment;

import au.edu.wehi.idsv.Defaults;
import com.google.common.io.Files;
import com.intel.gkl.smithwaterman.IntelSmithWaterman;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.UUID;

public class AlignerFactory {
	private static final Log log = Log.getInstance(AlignerFactory.class);
	private static final String SSW_JNI_JAR_LOCATION = "/libsswjni.so";
	private static final Aligner defaultAligner;
	private static boolean sswjniLoaded;
	private static final IntelSmithWaterman isw;
	static {
    	File tmpDir = new File(System.getProperty("java.io.tmpdir"));
    	if (!tmpDir.exists()) {
    		if (!tmpDir.mkdir()) {
				log.info("Created " + tmpDir);
			}
		}
    	log.debug("Loading Intel GKL library");
		IntelSmithWaterman initialisingIsw = new IntelSmithWaterman();
		if (Defaults.NO_LIBGKL || !initialisingIsw.load(tmpDir)) {
			initialisingIsw = null;
		} else {
			try {
				log.debug("Testing Intel GKL alignment");
				GKLAligner gklAligner = new GKLAligner(1, -4, 6, 1, initialisingIsw);
				gklAligner.align_smith_waterman(new byte[]{'A'}, new byte[]{'A'});
				log.info("Intel GKL library loading successful.");
			} catch (Exception e) {
				log.error(e, "Intel GKL call failure. Attempting to recover.");
				initialisingIsw = null;
			}
		}
		isw = initialisingIsw;
		if (isw == null) {
			if (!Defaults.NO_LIBGKL) {
				log.warn("Unable to use Intel GKL library for accelerated Smith-Waterman alignment");
			}
			if (!Defaults.NO_LIBSSW) {
				try {
					System.loadLibrary("sswjni");
					sswjniLoaded = true;
				} catch (UnsatisfiedLinkError e) {
					sswjniLoaded = false;
					log.debug(e, "Unable to load sswjni library");
				}
				try {
					if (!sswjniLoaded) {
						File tmp = new File(tmpDir, UUID.randomUUID().toString() + "-libsswjni.so");
						File destination = new File(tmpDir, "libsswjni.so");
						try {
							unpacksswjni(tmp);
							Files.move(tmp, destination);
							System.load(destination.getAbsolutePath());
							sswjniLoaded = true;
						} catch (IOException e) {
							log.error(e, "Unable to extract sswjni native library to " + destination.toString());
						}
					}
				} catch (UnsatisfiedLinkError e) {
					sswjniLoaded = false;
				}
			} else {
				sswjniLoaded = false;
			}
			if (sswjniLoaded) {
				try {
					log.debug("Testing ssw alignment");
					Aligner sswAligner = create(1, -4, -4, 6, 1);
					sswAligner.align_smith_waterman(new byte[]{'A'}, new byte[]{'A'});
					log.info("sswjni library loading successful.");
				} catch (Exception e) {
					log.error(e, "sswjni JNI call failure. Attempting to recover.");
					sswjniLoaded = false;
				}
			}
			if (!sswjniLoaded) {
				log.error("Unable to use GKL or sswjni libraries - realignment and inexact homology steps will be very slow. Please ensure Intel GKL and/or libsswjni for your OS and architecture can be found on java.library.path");
			}
		}
		// defaultAligner = create(2, -6, -1, 5, 3); // bowtie2 defaults
		defaultAligner = create(1, -4, -4, 6, 1); // bwa mem defaults
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
		if (isw != null) {
			if (mismatch != ambiguous) throw new IllegalArgumentException("GKL aligner does not support ambiguous scoring");
			return new GKLAligner(match, mismatch, gapOpen, gapExtend, isw);
		} else if (sswjniLoaded) {
			return new SswJniAligner(match, mismatch, ambiguous, gapOpen, gapExtend);
		} else {
			return new JAlignerAligner(match, mismatch, ambiguous, gapOpen, gapExtend);
		}
	}
	public static Aligner create() {
		return defaultAligner;
	}

	public static boolean isSswjniLoaded() { return sswjniLoaded; }
}
