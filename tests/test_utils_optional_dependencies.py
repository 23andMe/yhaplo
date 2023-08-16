import pytest

from yhaplo.utils.optional_dependencies import optional_import_error_message


def test_optional_import_error():
    with pytest.raises(ImportError, match="pip"):
        try:
            from uninstalled_module import some_function  # noqa F401
        except ImportError as error:
            error.msg = error.msg + optional_import_error_message(
                "package name",
                "do stuff",
                "optional dependency category",
            )
            raise error
