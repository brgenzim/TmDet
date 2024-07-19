<?php

namespace Unitmp\TmdetTest\Services;

use function PHPUnit\Framework\assertEquals;

class TmDetDsspRunner extends AbstractProcessRunner {

    const EXEC = '../build/src/tmdet';

    public array $chains = [];
    public array $dssps = [];

    protected function filterOutputLines(array $lines): array {

        $selectedLines = [];
        $beforeDebugLines = true;
        foreach($lines as $line) {
            if (str_starts_with($line, 'CHAIN ')) {
                $beforeDebugLines = false;
                list(, $chain, $residueCount) = explode(' ', $line);
                if ($residueCount !== "0" && !in_array($chain, $this->chains)) {
                    $this->chains[] = $chain;
                }
            }
            if ($beforeDebugLines && str_contains($line, ':')) {
                list($chain, $dsspString) = explode(':', $line);
                $selectedLines[] = trim($dsspString);
            }
        }

        $selectedLines = preg_grep('/^(\s*|.+at position.+)$/', $selectedLines, PREG_GREP_INVERT);
        if (count($this->chains) != count($selectedLines)) {
            printf("Command: %s\n", $this->commandLine);
            var_dump($this->chains);
            var_dump($selectedLines);
        }
        assertEquals(count($this->chains), count($selectedLines),
            'Chain counts do not match with number of dssp lines');
        $this->dssps = array_combine($this->chains, $selectedLines);

        return $selectedLines;
    }

    protected function parseOutputLine(string $line): string {
        // no action here
        return $line;
    }

    public static function createRunner(string $pdbFile): static {

        $pdbCode = static::getPdbCodeFromPath($pdbFile);
        $params = [
            '-i', "'$pdbFile'",
            '-x', "$pdbCode.xml"
        ];
        return new TmDetDsspRunner(static::EXEC, $params);
    }
}
