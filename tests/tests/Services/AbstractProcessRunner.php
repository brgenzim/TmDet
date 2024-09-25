<?php

namespace Unitmp\TmdetTest\Services;

use InvalidArgumentException;

abstract class AbstractProcessRunner {

    public readonly string $execPath;
    public readonly array $commandParams;
    public array $outputLines = [];
    protected string $commandLine;

    public function __construct(string $execPath, array $commandParams) {
        $this->execPath = $execPath;
        $this->commandParams = $commandParams;
    }

    /**
     * Returns subset of output lines. Any custom rule have to be implemented
     * in this function to get specific lines from the given array.
     */
    abstract protected function filterOutputLines(array $lines): array;

    /**
     * Parse a single line of the process output.
     */
    abstract protected function parseOutputLine(string $line): string;

    public function exec(): bool {
        $lines = null;
        $resultCode = null;

        $this->commandLine = "$this->execPath " . implode(' ', $this->commandParams);
        $this->commandLine .= ' 2>&1';
        exec($this->commandLine, $lines, $resultCode);

        if ($resultCode !== 0) {
            fprintf(STDERR, "Failed: %s\n", $this->commandLine);
            fprintf(STDERR, "Output:\n%s\n", implode("\n", $lines));
        }
        $lines = $this->filterOutputLines($lines);
        foreach ($lines as $line) {
                // just collect error message lines and jump to next line
                if ($resultCode !== 0) {
                $this->outputLines[] = $line;
                continue;
            }

            // collect output of successful run
            $this->outputLines[] = $this->parseOutputLine($line);
        }

        return $resultCode === 0;
    }

    public static function getPdbCodeFromPath(string $pdbPath): string {

        if (!str_ends_with($pdbPath, '.cif.gz') && !str_ends_with($pdbPath, '.ent.gz') && !str_ends_with($pdbPath, '.pdb.gz')) {
            throw new InvalidArgumentException("Unexpected file path: $pdbPath");
        }
        if (!preg_match('/(pdb)?([^\/]*?)\.(cif|ent|pdb)\.gz$/', $pdbPath, $matches)) {
            throw new InvalidArgumentException("Unexpected file path: $pdbPath");
        }

        return $matches[2];
    }
}
