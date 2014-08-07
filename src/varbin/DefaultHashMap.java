package varbin;

import java.util.HashMap;

public class DefaultHashMap<K,V> extends HashMap<K,V> {
	
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public V getDefault(K key, V defaultValue) {
        if (this.containsKey(key)) {
            return get(key);
        }
	        return defaultValue;
    }
}