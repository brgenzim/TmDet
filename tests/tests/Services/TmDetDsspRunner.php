<?php

namespace Unitmp\TmdetTest\Services;

class TmDetDsspRunner extends AbstractProcessRunner {

    const EXEC = '/home/A/csongor/dev/UniTmp/Tmdet/build/src/tmdet';

    public array $chains = [];
    public array $dssps = [];

    protected function filterOutputLines(array $lines): array {

        $selectedLines = [];
        $beforeDebugLines = true;
        foreach($lines as $line) {
            if (str_starts_with($line, 'CHAIN ')) {
                $beforeDebugLines = false;
                list(, $chain, ) = explode(' ', $line);
                $this->chains[] = $chain;
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

    public static function createRunner(string $pdbCode): static {
        $params = [
            '-i', "$pdbCode.cif.gz",
            '-x', "$pdbCode.xml"
        ];
        return new TmDetDsspRunner(static::EXEC, $params);
    }
}
