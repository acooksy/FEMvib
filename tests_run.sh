#!/bin/bash
        CYAN='\033[0;36m'
        GREEN='\033[0;32m'
        RED='\033[0;31m'
        NC='\033[0m' # No Color
        for test in tests/*; do
          if [ -e results/ ]; then
            rm -rf results/
          fi
          cp -pr $test ./results
          echo -e "${CYAN}++++++++  RUNNING TEST $test ++++++++${NC}"
          ./femvib_run.sh
          if diff -q results/eigenvalues.txt results/eigenvalues_std.txt; then
            echo -e "${GREEN}Test $test PASSED${NC}"
          else
            echo -e "${RED}Test $test FAILED${NC}"
            exit 1
          fi
        done
        echo -e "${GREEN}ALL TESTS PASSED${NC}"
