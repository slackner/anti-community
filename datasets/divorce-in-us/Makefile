all: metrics.txt

result.gml: result.graph
	../../src/gml "$^" > "$@"

metrics.txt: sbm.txt minmodularity.txt minmodularity-full.txt \
			 maxmodularity.txt laprop.txt antimodularity.txt
	@echo "" > metrics.txt
	for file in $^; do \
	    echo "[$${file%.txt}]" >> metrics.txt; \
	    ../../src/metrics --clusters="$$file" result.graph >> metrics.txt; \
	done

sbm.txt: result.gml
	../../third-party/estimate < result.gml > "$@"

minmodularity.txt: result.graph
	../../src/minmodularity --fast result.graph > minmodularity.txt

minmodularity-full.txt: result.graph
	../../src/minmodularity result.graph > minmodularity-full.txt

maxmodularity.txt: result.graph
	../../src/maxmodularity result.graph > maxmodularity.txt

laprop.txt: result.graph
	../../src/laprop result.graph > laprop.txt

antimodularity.txt: result.graph
	../../src/antimodularity result.graph > antimodularity.txt
