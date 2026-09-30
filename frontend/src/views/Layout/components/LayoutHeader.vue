<template>
  <header ref="header" class="site-header" @keydown.esc="onEscape" @focusout="onFocusOut">
    <nav class="container navigation" aria-label="Main navigation">
      <router-link to="/Home" class="brand" aria-label="SCAID home" @click="selectNavigation">
        <img src="@/assets/images/SCAID_nav.png" alt="SCAID" width="84" height="70" />
        <span class="brand-text">Single-Cell AutoImmune Disease database</span>
      </router-link>
      <button
        ref="menuButton"
        class="menu-toggle"
        type="button"
        :aria-expanded="menuOpen"
        aria-controls="scaid-navigation"
        @click="menuOpen = !menuOpen"
      >
        <svg aria-hidden="true" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.75" stroke-linecap="round">
          <path v-if="menuOpen" d="m6 6 12 12M18 6 6 18" />
          <path v-else d="M4 6h16M4 12h16M4 18h16" />
        </svg>
        <span>{{ menuOpen ? 'Close menu' : 'Menu' }}</span>
      </button>
      <div id="scaid-navigation" class="navigation-links" :class="{ 'is-open': menuOpen }">
        <ul class="navigation-list">
          <li v-for="item in navList" :key="item.link">
            <router-link
              :to="item.link"
              class="navigation-link"
              :aria-current="isCurrent(item.link) ? 'page' : undefined"
              @click="selectNavigation"
            >{{ item.name }}</router-link>
          </li>
        </ul>
      </div>
    </nav>
  </header>
</template>

<script setup>
import { nextTick, onBeforeUnmount, onMounted, ref, watch } from 'vue';
import { useRoute } from 'vue-router';

const route = useRoute();
const header = ref(null);
const menuButton = ref(null);
const menuOpen = ref(false);
const navList = [
  { name: 'Home', link: '/Home' },
  { name: 'Dataset', link: '/Dataset' },
  { name: 'Gene', link: '/Gene' },
  { name: 'KEGG', link: '/KEGG' },
  { name: 'Atlases', link: '/TeamInformation' },
  { name: 'Information', link: '/Information' },
];

const isCurrent = (path) => route.path.replace(/\/$/, '').toLowerCase() === path.toLowerCase();
const onEscape = (event) => {
  if (!menuOpen.value) return;
  event.preventDefault();
  menuOpen.value = false;
  menuButton.value?.focus();
};
const onFocusOut = (event) => {
  if (event.relatedTarget && !header.value?.contains(event.relatedTarget)) menuOpen.value = false;
};
const onPointerDown = (event) => {
  if (menuOpen.value && !header.value?.contains(event.target)) menuOpen.value = false;
};
const selectNavigation = async () => {
  const wasOpen = menuOpen.value;
  menuOpen.value = false;
  if (wasOpen) {
    await nextTick();
    document.getElementById('main-content')?.focus({ preventScroll: true });
  }
};
let desktopQuery;
const onBreakpointChange = (event) => {
  if (event.matches) menuOpen.value = false;
};

watch(() => route.fullPath, () => { menuOpen.value = false; });
onMounted(() => {
  document.addEventListener('pointerdown', onPointerDown);
  desktopQuery = window.matchMedia('(min-width: 992px)');
  desktopQuery.addEventListener('change', onBreakpointChange);
});
onBeforeUnmount(() => {
  document.removeEventListener('pointerdown', onPointerDown);
  desktopQuery?.removeEventListener('change', onBreakpointChange);
});
</script>

<style scoped lang="scss">
.site-header {
  position: sticky;
  top: 0;
  z-index: 100;
  width: 100%;
  background: #fff;
  border-top: 3px solid var(--scaid-accent, #8a503a);
  border-bottom: 1px solid var(--scaid-line, #e5ddd7);
}

.navigation {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 16px;
  min-height: 76px;
}

.brand {
  display: flex;
  align-items: center;
  flex-shrink: 0;
  gap: 12px;
  color: var(--scaid-muted, #6b5f5a);
  text-decoration: none;
  border-radius: 4px;
}

.brand img {
  display: block;
  width: 78px;
  height: 65px;
  object-fit: contain;
}

.brand-text {
  max-width: 290px;
  font-size: .875rem;
  line-height: 1.45;
}

.navigation-list {
  display: flex;
  align-items: center;
  gap: 4px;
  margin: 0;
  padding: 0;
  list-style: none;
}

.navigation-link {
  display: flex;
  align-items: center;
  min-height: 44px;
  padding: 10px 12px;
  border-radius: 4px;
  color: var(--scaid-ink, #2b2320);
  font-size: .9375rem;
  font-weight: 600;
  line-height: 1.5;
  white-space: nowrap;
  text-decoration: none;
  text-underline-offset: 8px;
  text-decoration-thickness: 2px;
  transition: color 160ms ease-out, background-color 160ms ease-out;
}

.navigation-link:hover,
.menu-toggle:hover {
  color: var(--scaid-accent, #8a503a);
  background: var(--scaid-soft, #f8f5f2);
}

.navigation-link[aria-current='page'] {
  color: var(--scaid-accent, #8a503a);
  text-decoration: underline;
}

.brand:focus-visible,
.navigation-link:focus-visible,
.menu-toggle:focus-visible {
  outline: 2px solid var(--scaid-accent, #8a503a);
  outline-offset: 3px;
}

.menu-toggle {
  display: none;
  align-items: center;
  justify-content: center;
  gap: 8px;
  min-height: 44px;
  padding: 9px 12px;
  border: 1px solid #c9bcb5;
  border-radius: 6px;
  color: var(--scaid-ink, #2b2320);
  background: #fff;
  font: inherit;
  font-size: .9375rem;
  font-weight: 600;
  line-height: 1.5;
  cursor: pointer;
}

@media (max-width: 1199px) {
  .brand-text { display: none; }
}

@media (max-width: 991px) {
  .navigation {
    flex-wrap: wrap;
    gap: 0 16px;
    min-height: 72px;
  }

  .brand img {
    width: 72px;
    height: 60px;
  }

  .menu-toggle { display: flex; }

  .navigation-links {
    display: none;
    flex-basis: 100%;
    max-height: calc(100dvh - 90px);
    overflow-y: auto;
    padding: 8px 0 12px;
    border-top: 1px solid var(--scaid-line, #e5ddd7);
  }

  .navigation-links.is-open { display: block; }

  .navigation-list {
    display: grid;
    gap: 2px;
  }

  .navigation-link {
    min-height: 46px;
    padding: 10px 12px;
    text-underline-offset: 4px;
  }

  .navigation-link[aria-current='page'] {
    background: var(--scaid-soft, #f8f5f2);
  }
}

@media (prefers-reduced-motion: reduce) {
  .navigation-link { transition: none; }
}
</style>
