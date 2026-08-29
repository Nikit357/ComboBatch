#!/bin/sh
# Container entry point.
#
# The plan called for `ENTRYPOINT ["combobatch"]`, but taking that literally would break
# every verification command that runs something else in the image - `pytest`, `Rscript`,
# `octave --version`, `python -c ...` - and the interactive pod, which overrides the
# command with `sleep infinity`. So anything that names a real executable runs as given,
# and everything else is handed to the CLI, which is what makes a mistyped subcommand an
# argparse error listing the real ones instead of `exec: nonsense: not found`.
#
#   docker run IMAGE                              -> combobatch --help   (the CMD)
#   docker run IMAGE combobatch selftest          -> combobatch selftest
#   docker run IMAGE --version                    -> combobatch --version
#   docker run IMAGE run --config c.yaml          -> combobatch run --config c.yaml
#   docker run IMAGE pytest tests/unit            -> pytest tests/unit
#
# Thread pinning is deliberately *not* done here. It belongs in the Python entry point,
# because a process launched over SSH into a running pod never passes through this file;
# the image's ENV and /etc/environment cover that case, and pin_threads() covers the rest.

set -e

if [ "$#" -eq 0 ]; then
    set -- combobatch --help
fi

case "$1" in
    -*) set -- combobatch "$@" ;;
    # No combobatch subcommand shares a name with a binary in this image, so "not an
    # executable" is a reliable test for "meant for the CLI".
    *) command -v "$1" >/dev/null 2>&1 || set -- combobatch "$@" ;;
esac

exec "$@"
