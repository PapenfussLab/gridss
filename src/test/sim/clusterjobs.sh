#!/bin/bash
# Author: Daniel Cameron
#
# check the server load on the WEHI computer servers
export TERM=xterm-256color

function check {
	echo $1 $(ssh -o ConnectTimeout=5 -o StrictHostKeyChecking=no -o ServerAliveCountMax=5 -o ServerAliveInterval=5 $1 ps -Fu cameron.d | grep time | wc -l )
}
for server in unix100.wehi.edu.au bionode01.wehi.edu.au bionode02.wehi.edu.au bionode03.wehi.edu.au bionode04.wehi.edu.au bionode05.wehi.edu.au bionode06.wehi.edu.au bionode07.wehi.edu.au bionode08.wehi.edu.au bionode09.wehi.edu.au bionode10.wehi.edu.au bionode11.wehi.edu.au bionode12.wehi.edu.au bionode13.wehi.edu.au bionode14.wehi.edu.au bionode15.wehi.edu.au bionode16.wehi.edu.au bionode17.wehi.edu.au bionode18.wehi.edu.au bionode19.wehi.edu.au bionode20.wehi.edu.au bionode21.wehi.edu.au bionode22.wehi.edu.au bionode23.wehi.edu.au bionode24.wehi.edu.au bionode25.wehi.edu.au bionode26.wehi.edu.au bionode27.wehi.edu.au bionode28.wehi.edu.au bionode29.wehi.edu.au bionode30.wehi.edu.au bionode31.wehi.edu.au bionode32.wehi.edu.au bionode33.wehi.edu.au ; do
	check $server &
done
wait

