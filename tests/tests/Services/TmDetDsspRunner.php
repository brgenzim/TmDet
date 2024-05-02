<?php

namespace Unitmp\TmdetTest\Services;

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
                list(, $chain, ) = explode(' ', $line);
                if (!in_array($chain, $this->chains)) {
                    $this->chains[] = $chain;
                }
            }
            if ($beforeDebugLines) {
                $selectedLines[] = $line;
            }
        }

        $selectedLines = preg_grep('/^\s*$/', $selectedLines, PREG_GREP_INVERT);
        $this->dssps = array_combine($this->chains, $selectedLines);

        return $selectedLines;
    }

    protected function parseOutputLine(string $line): string {
        // no action here
        return $line;
    }

    public static function createRunner(string $cifFile): static {

        $pdbCode = static::getPdbCodeFromZippedCifPath($cifFile);
        $params = [
            '-i', "'$cifFile'",
            '-x', "$pdbCode.xml"
        ];
        return new TmDetDsspRunner(static::EXEC, $params);
    }
}
