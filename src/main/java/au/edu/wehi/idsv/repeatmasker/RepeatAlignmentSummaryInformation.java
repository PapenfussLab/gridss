package au.edu.wehi.idsv.repeatmasker;

public class RepeatAlignmentSummaryInformation {
    private float percentageSubstituted;
    private float percentageDeleted;
    private float percentageInserted;
    private int basesInQueryPastMatch;
    private int matchStart;
    private int matchEnd;
    private int basesInRepeatPastMatch;

    public float getPercentageSubstituted() {
        return percentageSubstituted;
    }

    public float getPercentageDeleted() {
        return percentageDeleted;
    }

    public float getPercentageInserted() {
        return percentageInserted;
    }

    public int getBasesInQueryPastMatch() {
        return basesInQueryPastMatch;
    }

    public int getMatchStart() {
        return matchStart;
    }

    public int getMatchEnd() {
        return matchEnd;
    }

    public void setPercentageSubstituted(float percentageSubstituted) {
        this.percentageSubstituted = percentageSubstituted;
    }

    public void setPercentageDeleted(float percentageDeleted) {
        this.percentageDeleted = percentageDeleted;
    }

    public void setPercentageInserted(float percentageInserted) {
        this.percentageInserted = percentageInserted;
    }

    public void setBasesInQueryPastMatch(int basesInQueryPastMatch) {
        this.basesInQueryPastMatch = basesInQueryPastMatch;
    }

    public void setMatchStart(int matchStart) {
        this.matchStart = matchStart;
    }

    public void setMatchEnd(int matchEnd) {
        this.matchEnd = matchEnd;
    }

    public int getBasesInRepeatPastMatch() {
        return basesInRepeatPastMatch;
    }

    public void setBasesInRepeatPastMatch(int basesInRepeatPastMatch) {
        this.basesInRepeatPastMatch = basesInRepeatPastMatch;
    }
}
