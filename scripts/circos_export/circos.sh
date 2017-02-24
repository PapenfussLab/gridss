#!/bin/bash
#
# create a basic circos plot from a gridss BEDPE file
#
if [[ ! -f $1 ]] ; then
	echo "Error: input BEDPE not specified"
	exit 1
fi
if [[ ! -f $2 ]] ; then
	echo "Error: reference .dict not specified"
	exit 1
fi
BEDPE=$1
DICT=$2
SAMPLE=${1/.bedpe/}
mkdir -p $SAMPLE
bedpeToCircos.py < $1 > $SAMPLE/link.txt
dictToCircosKaryotype.py < $2 > $SAMPLE/karyotype.txt
cat > $SAMPLE/circos.conf << EOF
# circos.conf

<ideogram>
show_label     = yes
label_with_tag = yes
label_font     = light
label_radius   = 1.05r
label_center   = yes
label_size     = 12p
label_color    = grey
label_parallel = no
<spacing>
default = 0.001r
</spacing>
radius    = 0.9r
thickness = 20p
fill      = yes
</ideogram>


karyotype = karyotype.txt

chromosomes_display_default = yes

<links>
	<link>
		file = link.txt
		radius = 1r - 20p
		bezier_radius = 0.5r
		bezier_radius_purity = 0.1
	</link>
</links>

<image>
	<<include etc/image.conf>>
	auto_alpha_colors* = yes
	auto_alpha_steps* = 16
</image>

<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
EOF

( cd $SAMPLE && circos )

