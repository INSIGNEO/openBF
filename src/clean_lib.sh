#!/bin/bash
find . -type f -not -name '*test*' | xargs rm -rf
rm -rf doc/
