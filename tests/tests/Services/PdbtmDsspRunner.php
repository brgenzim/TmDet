<?php

namespace Unitmp\TmdetTest\Services;

class PdbtmDsspRunner extends AbstractProcessRunner {

    const EXEC = '/home/tusi/works/pdbtm_3.0/TmdetUtils/bin/dssp';
    const HEADER_LINE_PREFIX = '  #  RESIDUE AA STRUCTURE BP1 BP2  ACC';
    const CHAIN_DELIMITER = '      !';
    const CHAIN_COLUMN = 1;
    const STRUCTURE_COLUMN = 5;

    public array $chains = [];
    public array $dssps = [];
    public string $dsspCacheFile = '';

    protected function filterOutputLines(array $lines): array {

        // first 50 lines should be enough to get dssp file path
        $dsspMatches = preg_grep('/^\/.+[0-9]+.dssp\s*$/', array_slice($lines, 0, 50));
        if (!empty($dsspMatches) && file_exists($dsspFile = array_shift($dsspMatches))) {
            $this->dsspCacheFile = $dsspFile;

            return explode("\n", file_get_contents($this->dsspCacheFile));
        }

        $selectedLines = [];
        $recordLines = false;
        foreach($lines as $line) {
            if (str_starts_with($line, static::HEADER_LINE_PREFIX)) {
                $recordLines = true;
                continue;
            }
            if ($recordLines) {
                $selectedLines[] = $line;
            }
        }

        return $selectedLines;
    }

    protected function parseOutputLine(string $line): string {

        if (!empty($this->dsspCacheFile)) {
            return $this->parseDsspLine($line);
        }

        if (str_contains($line, static::CHAIN_DELIMITER)) {
            return $line;
        }

        // TODO: collect dssp strings when there is no cache file

        return $line;
    }

    protected function parseDsspLine(string $line): string {
        // parse a chain line
        if (preg_match('/^(\w+)_(\w+)/', $line, $matches)) {
            $chain = $matches[2];
            if (!in_array($chain, $this->chains)) {
                $this->chains[] = $chain;
            }
        } elseif (count(explode(' ', $line)) == 7) {
            $column = $line[static::CHAIN_COLUMN - 1];
            if (!array_key_exists($column, $this->dssps)) {
                $this->dssps[$column] = '';
            }
            $secondaryStructure = $line[static::STRUCTURE_COLUMN - 1];
            if ($secondaryStructure == ' ') {
                $secondaryStructure = '-';
            }
            $this->dssps[$column] = $this->dssps[$column]
                . $secondaryStructure;
        }
        return $line;
    }

    public static function createRunner(string $cifFile): static {

        $params = [
            '-i', "'$cifFile'",
            '--output-format', 'dssp'
        ];
        return new PdbtmDsspRunner(static::EXEC, $params);
    }
}
