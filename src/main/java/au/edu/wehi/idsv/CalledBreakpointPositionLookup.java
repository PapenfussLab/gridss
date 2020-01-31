package au.edu.wehi.idsv;

import java.util.HashMap;

public class CalledBreakpointPositionLookup {
    public static class NominalPosition {
        public final BreakpointSummary nominalPosition;
        public final String insertedSequenced;
        public final String homologySequence;
        public final boolean isExact;

        public NominalPosition(BreakpointSummary nominalPosition, String insertedSequenced, String homologySequence, boolean isExact) {
            this.nominalPosition = nominalPosition;
            this.insertedSequenced = insertedSequenced;
            this.homologySequence = homologySequence;
            this.isExact = isExact;
        }

        public NominalPosition remoteBreakpoint() {
            return new NominalPosition(
                    nominalPosition.remoteBreakpoint(),
                    nominalPosition.direction != nominalPosition.direction2 ? insertedSequenced : htsjdk.samtools.util.SequenceUtil.reverseComplement(insertedSequenced),
                    nominalPosition.direction != nominalPosition.direction2 ? homologySequence : htsjdk.samtools.util.SequenceUtil.reverseComplement(homologySequence),
                    isExact
            );
        }
    }
    private HashMap<String, NominalPosition> lookup = new HashMap<>();
    public void addLower(String eventId, NominalPosition position) {
        lookup.put(eventId, position.remoteBreakpoint());
    }
    public NominalPosition removeUpper(String eventId) {
        return lookup.remove(eventId);
    }
}
