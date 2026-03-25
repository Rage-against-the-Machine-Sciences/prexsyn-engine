#!/bin/bash

rm -rf prexsyn_engine-stubs
stubgen -p prexsyn_engine -o typings
mv typings/prexsyn_engine prexsyn_engine-stubs
