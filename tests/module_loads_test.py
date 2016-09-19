def test_module_loads():
    loaded = False
    try:
        import yhaplo
        loaded = True
    except ImportError:
        pass

    assert loaded
