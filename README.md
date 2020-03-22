# QDM

## Building from Source

Build the binary using docker:

```bash
git clone https://github.com/calebcase/qdm.git
cd qdm
docker build -t qdm .
```

Optionally, extract the binary from the docker image:

```bash
docker create -it --name qdm qdm
docker cp qdm:/usr/local/bin/qdm qdm
docker rm -f qdm
```

## Running

Run the binary from within docker:

```bash
docker run -it --rm -v $(pwd):/data qdm
```

Run locally (if you have extracted the binary):

```bash
./qdm
```
