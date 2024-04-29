<?php

namespace Unitmp\TmdetTest\Services;

class DsspRunner extends AbstractProcessRunner {

    const EXEC = 'dssp';
    const HEADER_LINE_PREFIX = '  #  RESIDUE AA STRUCTURE BP1 BP2  ACC';
    const CHAIN_DELIMITER = '      !';
    const CHAIN_COLUMN = 12;
    const STRUCTURE_COLUMN = 17;

    public array $chains = [];
    public array $dssps = [];

    protected function filterOutputLines(array $lines): array {

        $selectedLines = [];
        $recordLines = false;
        foreach($lines as $line) {
            if (str_starts_with($line, static::HEADER_LINE_PREFIX)) {
                $recordLines = true;
            }
            if ($recordLines) {
                $selectedLines[] = $line;
            }
        }

        return $selectedLines;
    }

    protected function parseOutputLine(string $line): string {

        if (str_contains($line, static::CHAIN_DELIMITER)) {
            return $line;
        }

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

        return $line;
    }

    public static function createRunner(string $pdbCode): static {
        $params = [
            '-i', "$pdbCode.cif.gz",
            '--output-format', 'dssp'
        ];
        return new DsspRunner(static::EXEC, $params);
    }
}
