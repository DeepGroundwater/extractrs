from pipelines.channel_attrs import paths


def test_inputs_exist():
    assert paths.MERIT_RIV.exists(), "MERIT CONUS flowlines missing"
    assert paths.GLOBAL_NC.exists(), "global attributes nc missing"


def test_staging_dirs():
    assert paths.RAW.is_dir() and paths.DERIVED.is_dir()
