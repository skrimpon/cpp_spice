#!/bin/bash
find . -name '*.h' -o -name '*.cpp' | xargs wc -l
