#!/bin/bash
# Based on code by @gittywithexcitement :
# https://github.com/commercialhaskell/stack/issues/5228#issuecomment-629738842

stack haddock

LOCAL_DOC_ROOT=$(stack path --local-doc-root)
SNAPSHOT_DOC=$(stack path --snapshot-doc-root)
DATABASE="$(stack path --local-hoogle-root)/database.foo"

stack exec -- hoogle generate --local=$LOCAL_DOC_ROOT --local=$SNAPSHOT_DOC --database=$DATABASE

hoogle server --local --database=$DATABASE > /dev/null &

