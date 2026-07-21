#!/bin/bash

echo "==== User Permissions Report ===="
for user in $(cut -d: -f1 /etc/passwd); do
    home_dir=$(eval echo ~$user)
    echo "User: $user"
    
    # Home directory permissions
    if [ -d "$home_dir" ]; then
        echo "Home dir: $home_dir"
        ls -ld "$home_dir"
    fi
    
    # Group memberships
    echo -n "Groups: "
    groups $user

    # Sudo privileges
    echo "Sudo privileges:"
    sudo -lU $user 2>/dev/null | grep -v "may run the following"

    echo "----------------------------------"
done
