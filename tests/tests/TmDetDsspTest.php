<?php

namespace Unitmp\TmdetTest;

use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\Attributes\Group;
use PHPUnit\Framework\TestCase;
use Unitmp\TmdetTest\Services\DsspRunner;
use Unitmp\TmdetTest\Services\TmDetDsspRunner;

class TmDetDsspTest extends TestCase {

    const PDB_ZFS_DIR = '/zfs/databases/UniTmp/PDB/data/structures/divided/mmCIF';

    #[Group('dssp')]
    public function test_1afo_dssp() {
        // Arrange
        $tmdetRunner = TmDetDsspRunner::createRunner('1AFO');
        $dsspRunner = DsspRunner::createRunner('1AFO');
        // Act
        $tmdetSuccess = $tmdetRunner->exec();
        $dsspSuccess = $dsspRunner->exec();
        // Assert
        $this->assertTrue($tmdetSuccess);
        $this->assertTrue($dsspSuccess);
    }

    #[Group('dssp')]
    public function test_7rh7_dssp() {
        // Arrange
        $tmdetRunner = TmDetDsspRunner::createRunner('7RH7');
        $dsspRunner = DsspRunner::createRunner('7RH7');
        // Act
        $dsspSuccess = $dsspRunner->exec();
        $tmdetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($tmdetSuccess);
        $this->assertTrue($dsspSuccess);
        $this->assertEquals($dsspRunner->dssps['E'], $tmdetRunner->dssps['E']);
        $this->assertEquals($dsspRunner->dssps, $tmdetRunner->dssps);
    }

    #[Group('dssp')]
    #[DataProvider('getCifPaths')]
    public function test_whole_archive(string $cifPath) {

        // Pre-check
        if (static::isStructureUnsupported($cifPath)) {
            $this->markTestSkipped('DNA/RNA structures are unsupported');
            return;
        }

        // Arrange
        $tmdetRunner = TmDetDsspRunner::createRunner($cifPath);
        $dsspRunner = DsspRunner::createRunner($cifPath);
        // Act
        $dsspSuccess = $dsspRunner->exec();
        $tmdetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($tmdetSuccess);
        $this->assertTrue($dsspSuccess);
        $this->assertEquals($dsspRunner->dssps, $tmdetRunner->dssps);
    }

    //
    // Util functions
    //

    /**
     * Data provider for test_whole_archive
     */
    public static function getCifPaths(): array {

        $dirs = scandir(static::PDB_ZFS_DIR);
        if (!$dirs || count($dirs) == 2) {
            return [];
        }

        $files = [];
        foreach (array_slice($dirs, 2) as $directory) {
            $directory = static::PDB_ZFS_DIR . "/$directory";
            $items = scandir($directory);
            $items = preg_grep('/^.+cif.gz$/', $items);
            foreach($items as $item) {
                $files[ $item ] = [ "$directory/$item" ];
            }
        }
        // TODO: remove slice later
        return array_slice($files, 0, 200);
    }

    public static function isStructureUnsupported(string $cifPath): bool {
        exec("zgrep -P 'DNA|RNA' '$cifPath'", $lines, $resultCode);
        return $resultCode === 0;
    }
}
