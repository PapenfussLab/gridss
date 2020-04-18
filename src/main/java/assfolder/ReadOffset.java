package assfolder;

public class ReadOffset {
    public final Read read;
    public final int offset;
    public ReadOffset(Read read, int offset) {
        this.read = read;
        this.offset = offset;
    }

    @Override
    public boolean equals(Object obj) {
        ReadOffset ro = (ReadOffset)obj;
        return ro.read == read && ro.offset == offset;
    }

    @Override
    public int hashCode() {
        return offset + read.hashCode();
    }
}
