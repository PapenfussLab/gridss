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
 * Metrics about the distribution of extracted structural variant reads 
 *
 * @author Daniel Cameron
 */
public class StructuralVariantReadMetrics extends MultilevelMetrics {
	public long STRUCTURAL_VARIANT_READS;
	public long STRUCTURAL_VARIANT_READ_PAIRS;
    public long INDEL_READS;
    public long SPLIT_READS;
    public long SOFT_CLIPPED_READS;
    public long UNMAPPED_READS;
    public long DISCORDANT_READ_PAIRS;
    public long UNMAPPED_MATE_READS;
    public long STRUCTURAL_VARIANT_READ_ALIGNMENTS;
    public long INDEL_READ_ALIGNMENTS;
    public long SPLIT_READ_ALIGNMENTS;
    public long SOFT_CLIPPED_READ_ALIGNMENTS;
    /**
     * Read alignments in which the mate is discordant. libraries containing singlely mapped reads,
     * this value is twice DISCORDANT_READ_PAIRS.
     */
    public long DISCORDANT_READ_PAIR_ALIGNMENTS;
    public long UNMAPPED_MATE_READ_ALIGNMENTS;
}
