package au.edu.wehi.idsv.sim;

public interface RandomSequenceGenerator {
    byte[] getBases(int length);
    default Integer lastStartPosition() { return null; }
    default Integer lastEndPosition() { return null; }
}
