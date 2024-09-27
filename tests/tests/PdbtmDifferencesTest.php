<?php

namespace Unitmp\TmdetTest;

use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\Attributes\Group;
use PHPUnit\Framework\TestCase;
use Unitmp\TmdetTest\Constants\FileSystem;
use Unitmp\TmdetTest\Services\Tmdet30Runner;

class PdbtmDifferencesTest extends TestCase {

    const PDB_ENT_DIR = FileSystem::PDB_ENT_ZFS_DIR;
    const MAX_ALLOWED_LEVENSHTEIN_DISTANCE_IN_PERCENT = 40;

    #[Group('tmdet')]
    #[DataProvider('getEntPathsOfTmps')]
    public function test_pdbtm_manual_updates(string $entPath) {
        // Arrange
        $tmdetRunner = Tmdet30Runner::createRunner($entPath);
        // Act
        $tmdetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($tmdetSuccess);
        $this->assertRegionArraysEqual($tmdetRunner->oldData, $tmdetRunner->newData);
    }

    //
    // Util functions
    //
    public function assertRegionArraysEqual(array $expected, array $actual) {
        $this->assertEquals(array_keys($expected), array_keys($actual));
        $ok = true;
        foreach ($expected as $key => $expectedRegionString) {
            $actualRegionString = $actual[$key];
            if (strlen($expectedRegionString) != strlen($actualRegionString)) {
                break;
            }
            // max 2% difference allowed
            $maxDistance = round(strlen($expectedRegionString) / 100 * static::MAX_ALLOWED_LEVENSHTEIN_DISTANCE_IN_PERCENT, 2);
            $actualDistance = levenshtein($expectedRegionString, $actualRegionString);
            $actualPercent = round($actualDistance / strlen($actualRegionString) * 100, 2);
            if ($actualDistance > $maxDistance) {
                $ok = false;
                break;
            }
        }

        if (!isset($expectedRegionString)) {
            $this->assertTrue(true);
            return;
        }

        $this->assertEquals(strlen($expectedRegionString), strlen($actualRegionString),
            'Region compare failed at chain ' . $key . ' - region string lengths differ'
            . "\nexpected: '$expectedRegionString'"
            . "\nactual:   '$actualRegionString'");

        $this->assertTrue($ok, 'Region compare failed at chain ' . $key
            . "\nexpected: '$expectedRegionString'"
            . "\nactual:   '$actualRegionString'"
            . "\nmax expected levenshtein distance (%): "
            .  static::MAX_ALLOWED_LEVENSHTEIN_DISTANCE_IN_PERCENT
            . ", actual distance (%): $actualPercent"
            . "\nmax expected levenshtein distance: $maxDistance, actual distance: $actualDistance");
    }

    public static function getEntPathsOfTmps(): array {
        $paths = [];
        foreach (explode("\n", trim(file_get_contents('../data/tm-codes.txt'))) as $code) {
            $subDir = '/' . $code[1] . $code[2] . '/';
            $fileName = "pdb$code.ent.gz";
            $fullPath = static::PDB_ENT_DIR . $subDir . $fileName;
            $paths[ $fileName ] = [ $fullPath ];
        }
        return $paths;
        // return array_slice($paths, 0, 5);
    }

}
