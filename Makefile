SHELL := bash

all: qdm

.PHONY: image
image:
	docker build -t qdm .

qdm: image
	docker create -it --name qdm qdm
	docker cp qdm:/usr/local/bin/qdm qdm
	docker rm -f qdm
