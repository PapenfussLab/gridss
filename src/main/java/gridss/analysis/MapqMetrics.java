/*
 * The MIT License
 *
 * Copyright (c) 2016 Daniel Cameron
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package gridss.analysis;

import picard.metrics.MultilevelMetrics;

/**
 * Metrics about the mapping quality distribution of a paired-end library, created by the
 * CollectMapqMetrics program and usually written to a file with the extension
 * ".mapq_metrics".  In addition the insert size distribution is plotted to
 * a file with the extension ".mapq_Histogram.pdf".
 *
 * @author Daniel Cameron
 */
public class MapqMetrics extends MultilevelMetrics {
    
    public int MAX_MAPQ;
    
    public int MIN_MAPQ;
    
    public long MAPPED_READS;
    
    public long UNKNOWN_MAPQ;
    
    public long ZERO_MAPQ;
}
