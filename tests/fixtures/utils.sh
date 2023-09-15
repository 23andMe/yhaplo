function echo_run {
    # Echo and run a command.
    # Usage: echo_run "command [options] [arguments]"

    command=$1
    echo -e "\n$ ${command}"
    eval ${command}
}
