<?php

namespace Unitmp\TmdetTest\Services;

class PdbtmDsspRunner extends AbstractProcessRunner {

    const EXEC = 'export PDB_CACHE_PATH=/tmp; /home/tusi/works/pdbtm_3.0/TmdetUtils/bin/dssp';
    //const EXEC = 'export PDB_CACHE_PATH=/tmp; export LD_LIBRARY_PATH=/home/A/csongor/dev/PdbLib.C__stale/lib; /home/A/csongor/dev/TmdetUtils/bin/dssp';
    const DUMP_HEADER_LINE_PREFIX = 'ShowProt:';
    const CHAIN_COLUMN = 1;
    const STRUCTURE_COLUMN = 5;

    public array $chains = [];
    public array $dssps = [];
    public string $dsspCacheFile = '';

    protected function filterOutputLines(array $lines): array {

        $lines = preg_grep('/^PDB_WARNING:/', $lines, PREG_GREP_INVERT);
        // first 50 lines should be enough to get dssp file path
        $dsspMatches = preg_grep('/^\/.+[0-9]+.dssp\s*$/', array_slice($lines, 0, 50));
        if (!empty($dsspMatches) && file_exists($dsspFile = array_shift($dsspMatches))) {
            $this->dsspCacheFile = $dsspFile;

            return explode("\n", file_get_contents($this->dsspCacheFile));
        }

        $selectedLines = [];
        $recordLines = true;
        foreach($lines as $line) {

            if (preg_match('/\s+CHAIN\s+(\S+)\s+(\d+)/', $line, $matches)
                && $matches[2] !== "0") {

                $this->chains[] = $matches[1];
            }

            if (str_starts_with($line, static::DUMP_HEADER_LINE_PREFIX)) {
                $recordLines = false;
                continue;
            }
            if ($recordLines && preg_match('/^[-BEGHIST]+$/', $line, $matches)) {
                $selectedLines[] = $line;
            }
        }
        if (count($this->chains) != count($selectedLines)) {
            printf("Command: %s\n", $this->commandLine);
            var_dump($this->chains);
            var_dump($selectedLines);
        }
        $this->dssps = array_combine($this->chains, $selectedLines);

        return $selectedLines;
    }

    protected function parseOutputLine(string $line): string {

        if (!empty($this->dsspCacheFile)) {
            return $this->parseDsspLine($line);
        }

        return $line;
    }

    protected function parseDsspLine(string $line): string {
        // parse a chain line
        if (preg_match('/^(\w+)_(\w+)/', $line, $matches)) {
            $chain = $matches[2];
            if (!in_array($chain, $this->chains)) {
                $this->chains[] = $chain;
            }
        } elseif (count(($columns = explode(' ', $line))) == 7) {
            $column = $columns[static::CHAIN_COLUMN - 1];
            if (!array_key_exists($column, $this->dssps)) {
                $this->dssps[$column] = '';
            }
            $secondaryStructure = $columns[static::STRUCTURE_COLUMN - 1];
            if ($secondaryStructure == ' ') {
                $secondaryStructure = '-';
            }
            $this->dssps[$column] = $this->dssps[$column]
                . $secondaryStructure;
        }
        return $line;
    }

    public static function createRunner(string $entFile): static {

        $params = [
            '-i', "'$entFile'",
            '--output-format', 'dssp'
        ];
        return new PdbtmDsspRunner(static::EXEC, $params);
    }
}
