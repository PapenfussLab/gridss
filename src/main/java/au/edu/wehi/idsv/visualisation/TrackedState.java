package au.edu.wehi.idsv.visualisation;

import java.util.Collection;

public interface TrackedState {
    String[] trackedNames();
    Object[] trackedState();
    Collection<TrackedState> trackedObjects();
}
