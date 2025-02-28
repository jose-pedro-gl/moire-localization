function hash = hash_string(string)
    string = convertStringsToChars(string);
    md = java.security.MessageDigest.getInstance("SHA-256");
    hash = sprintf("%02X", typecast(md.digest(uint8(string)), "uint8")');
end
