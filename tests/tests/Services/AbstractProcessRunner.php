<?php

namespace Unitmp\TmdetTest\Services;

abstract class AbstractProcessRunner {

    public readonly string $execPath;
    public readonly array $commandParams;
    public array $outputLines = [];

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

        $command = "$this->execPath " . implode(' ', $this->commandParams);
        exec($command, $lines, $resultCode);

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
}
